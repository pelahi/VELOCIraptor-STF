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

/// Once data is loaded using field names in HDF5 file, update the names
/// to be unique to the selection loaded by VR
inline void UpdateExtraFieldNames(Options &opt)
{
    for (auto i=0;i<opt.gas_internalprop_names.size();i++) {
        opt.gas_internalprop_names[i]+= to_string(opt.gas_internalprop_index[i]);
    }
    for (auto i=0;i<opt.gas_chem_names.size();i++) {
        opt.gas_chem_names[i]+= to_string(opt.gas_chem_index[i]);
    }
    for (auto i=0;i<opt.gas_chemproduction_names.size();i++) {
        opt.gas_chemproduction_names[i]+= to_string(opt.gas_chemproduction_index[i]);
    }
    for (auto i=0;i<opt.star_internalprop_names.size();i++) {
        opt.star_internalprop_names[i]+= to_string(opt.star_internalprop_index[i]);
    }
    for (auto i=0;i<opt.star_chem_names.size();i++) {
        opt.star_chem_names[i]+= to_string(opt.star_chem_index[i]);
    }
    for (auto i=0;i<opt.star_chemproduction_names.size();i++) {
        opt.star_chemproduction_names[i]+= to_string(opt.star_chemproduction_index[i]);
    }
    for (auto i=0;i<opt.bh_internalprop_names.size();i++) {
        opt.bh_internalprop_names[i]+= to_string(opt.bh_internalprop_index[i]);
    }
    for (auto i=0;i<opt.bh_chem_names.size();i++) {
        opt.bh_chem_names[i]+= to_string(opt.bh_chem_index[i]);
    }
    for (auto i=0;i<opt.bh_chemproduction_names.size();i++) {
        opt.bh_chemproduction_names[i]+= to_string(opt.bh_chemproduction_index[i]);
    }
    for (auto i=0;i<opt.extra_dm_internalprop_names.size();i++) {
        opt.extra_dm_internalprop_names[i]+= to_string(opt.extra_dm_internalprop_index[i]);
    }
    for (auto i=0;i<opt.gas_internalprop_names_aperture.size();i++) {
        opt.gas_internalprop_names_aperture[i]+= to_string(opt.gas_internalprop_index_aperture[i]);
    }
    for (auto i=0;i<opt.gas_chem_names_aperture.size();i++) {
        opt.gas_chem_names_aperture[i]+= to_string(opt.gas_chem_index_aperture[i]);
    }
    for (auto i=0;i<opt.gas_chemproduction_names_aperture.size();i++) {
        opt.gas_chemproduction_names_aperture[i]+= to_string(opt.gas_chemproduction_index_aperture[i]);
    }
    for (auto i=0;i<opt.star_internalprop_names_aperture.size();i++) {
        opt.star_internalprop_names_aperture[i]+= to_string(opt.star_internalprop_index_aperture[i]);
    }
    for (auto i=0;i<opt.star_chem_names_aperture.size();i++) {
        opt.star_chem_names_aperture[i]+= to_string(opt.star_chem_index_aperture[i]);
    }
    for (auto i=0;i<opt.star_chemproduction_names_aperture.size();i++) {
        opt.star_chemproduction_names_aperture[i]+= to_string(opt.star_chemproduction_index_aperture[i]);
    }
    for (auto i=0;i<opt.bh_internalprop_names_aperture.size();i++) {
        opt.bh_internalprop_names_aperture[i]+= to_string(opt.bh_internalprop_index_aperture[i]);
    }
    for (auto i=0;i<opt.bh_chem_names_aperture.size();i++) {
        opt.bh_chem_names_aperture[i]+= to_string(opt.bh_chem_index_aperture[i]);
    }
    for (auto i=0;i<opt.bh_chemproduction_names_aperture.size();i++) {
        opt.bh_chemproduction_names_aperture[i]+= to_string(opt.bh_chemproduction_index_aperture[i]);
    }
    for (auto i=0;i<opt.extra_dm_internalprop_names_aperture.size();i++) {
        opt.extra_dm_internalprop_names_aperture[i]+= to_string(opt.extra_dm_internalprop_index_aperture[i]);
    }
}

inline void SetUniqueInputNames(Options & opt){
#ifdef GASON
        set<string> unique_gas_internalprop_names, unique_gas_chem_names, unique_gas_chemproduction_names;
        for (auto iextra=0;iextra<opt.gas_internalprop_names.size();iextra++) {
            string s = opt.gas_internalprop_names[iextra]+to_string(opt.gas_internalprop_index[iextra]);
            if (unique_gas_internalprop_names.count(s)==0) {
                unique_gas_internalprop_names.insert(s);
                opt.gas_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.gas_internalprop_names_aperture.size();iextra++) {
            string s = opt.gas_internalprop_names_aperture[iextra]+to_string(opt.gas_internalprop_index_aperture[iextra]);
            if (unique_gas_internalprop_names.count(s)==0) {
                unique_gas_internalprop_names.insert(s);
                opt.gas_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++) {
            string s = opt.gas_chem_names[iextra]+to_string(opt.gas_chem_index[iextra]);
            if (unique_gas_chem_names.count(s)==0) {
                unique_gas_chem_names.insert(s);
                opt.gas_chem_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.gas_chem_names_aperture.size();iextra++) {
            string s = opt.gas_chem_names_aperture[iextra]+to_string(opt.gas_chem_index_aperture[iextra]);
            if (unique_gas_chem_names.count(s)==0) {
                unique_gas_chem_names.insert(s);
                opt.gas_chem_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.gas_chemproduction_names.size();iextra++) {
            string s = opt.gas_chemproduction_names[iextra]+to_string(opt.gas_chemproduction_index[iextra]);
            if (unique_gas_chemproduction_names.count(s)==0) {
                unique_gas_chemproduction_names.insert(s);
                opt.gas_chemproduction_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.gas_chemproduction_names_aperture.size();iextra++) {
            string s = opt.gas_chemproduction_names_aperture[iextra]+to_string(opt.gas_chemproduction_index_aperture[iextra]);
            if (unique_gas_chemproduction_names.count(s)==0) {
                unique_gas_chemproduction_names.insert(s);
                opt.gas_chemproduction_unique_input_indexlist.push_back(iextra);
            }
        }
        opt.gas_internalprop_unique_input_names.assign(unique_gas_internalprop_names.begin(), unique_gas_internalprop_names.end());
        opt.gas_chem_unique_input_names.assign(unique_gas_chem_names.begin(), unique_gas_chem_names.end());
        opt.gas_chemproduction_unique_input_names.assign(unique_gas_chemproduction_names.begin(), unique_gas_chemproduction_names.end());
#endif
#ifdef STARON
        set<string> unique_star_internalprop_names, unique_star_chem_names, unique_star_chemproduction_names;
        for (auto iextra=0;iextra<opt.star_internalprop_names.size();iextra++) {
            string s = opt.star_internalprop_names[iextra]+to_string(opt.star_internalprop_index[iextra]);
            if (unique_star_internalprop_names.count(s)==0) {
                unique_star_internalprop_names.insert(s);
                opt.star_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.star_internalprop_names_aperture.size();iextra++) {
            string s = opt.star_internalprop_names_aperture[iextra]+to_string(opt.star_internalprop_index_aperture[iextra]);
            if (unique_star_internalprop_names.count(s)==0) {
                unique_star_internalprop_names.insert(s);
                opt.star_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.star_chem_names.size();iextra++) {
            string s = opt.star_chem_names[iextra]+to_string(opt.star_chem_index[iextra]);
            if (unique_star_chem_names.count(s)==0) {
                unique_star_chem_names.insert(s);
                opt.star_chem_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.star_chem_names_aperture.size();iextra++) {
            string s = opt.star_chem_names_aperture[iextra]+to_string(opt.star_chem_index_aperture[iextra]);
            if (unique_star_chem_names.count(s)==0) {
                unique_star_chem_names.insert(s);
                opt.star_chem_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.star_chemproduction_names.size();iextra++) {
            string s = opt.star_chemproduction_names[iextra]+to_string(opt.star_chemproduction_index[iextra]);
            if (unique_star_chemproduction_names.count(s)==0) {
                unique_star_chemproduction_names.insert(s);
                opt.star_chemproduction_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.star_chemproduction_names_aperture.size();iextra++) {
            string s = opt.star_chemproduction_names_aperture[iextra]+to_string(opt.star_chemproduction_index_aperture[iextra]);
            if (unique_star_chemproduction_names.count(s)==0) {
                unique_star_chemproduction_names.insert(s);
                opt.star_chemproduction_unique_input_indexlist.push_back(iextra);
            }
        }
        opt.star_internalprop_unique_input_names.assign(unique_star_internalprop_names.begin(), unique_star_internalprop_names.end());
        opt.star_chem_unique_input_names.assign(unique_star_chem_names.begin(), unique_star_chem_names.end());
        opt.star_chemproduction_unique_input_names.assign(unique_star_chemproduction_names.begin(), unique_star_chemproduction_names.end());
#endif
#ifdef BHON
        set<string> unique_bh_internalprop_names, unique_bh_chem_names, unique_bh_chemproduction_names;
        for (auto iextra=0;iextra<opt.bh_internalprop_names.size();iextra++) {
            string s = opt.bh_internalprop_names[iextra]+to_string(opt.bh_internalprop_index[iextra]);
            if (unique_bh_internalprop_names.count(s)==0) {
                unique_bh_internalprop_names.insert(s);
                opt.bh_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.bh_internalprop_names_aperture.size();iextra++) {
            string s = opt.bh_internalprop_names_aperture[iextra]+to_string(opt.bh_internalprop_index_aperture[iextra]);
            if (unique_bh_internalprop_names.count(s)==0) {
                unique_bh_internalprop_names.insert(s);
                opt.bh_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.bh_chem_names.size();iextra++) {
            string s = opt.bh_chem_names[iextra]+to_string(opt.bh_chem_index[iextra]);
            if (unique_bh_chem_names.count(s)==0) {
                unique_bh_chem_names.insert(s);
                opt.bh_chem_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.bh_chem_names_aperture.size();iextra++) {
            string s = opt.bh_chem_names_aperture[iextra]+to_string(opt.bh_chem_index_aperture[iextra]);
            if (unique_bh_chem_names.count(s)==0) {
                unique_bh_chem_names.insert(s);
                opt.bh_chem_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.bh_chemproduction_names.size();iextra++) {
            string s = opt.bh_chemproduction_names[iextra]+to_string(opt.bh_chemproduction_index[iextra]);
            if (unique_bh_chemproduction_names.count(s)==0) {
                unique_bh_chemproduction_names.insert(s);
                opt.bh_chemproduction_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.bh_chemproduction_names_aperture.size();iextra++) {
            string s = opt.bh_chemproduction_names_aperture[iextra]+to_string(opt.bh_chemproduction_index_aperture[iextra]);
            if (unique_bh_chemproduction_names.count(s)==0) {
                unique_bh_chemproduction_names.insert(s);
                opt.bh_chemproduction_unique_input_indexlist.push_back(iextra);
            }
        }
        opt.bh_internalprop_unique_input_names.assign(unique_bh_internalprop_names.begin(), unique_bh_internalprop_names.end());
        opt.bh_chem_unique_input_names.assign(unique_bh_chem_names.begin(), unique_bh_chem_names.end());
        opt.bh_chemproduction_unique_input_names.assign(unique_bh_chemproduction_names.begin(), unique_bh_chemproduction_names.end());
#endif
#ifdef EXTRADMON
        set<string> unique_extra_dm_internalprop_names, unique_extra_dm_chem_names, unique_extra_dm_chemproduction_names;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_names.size();iextra++) {
            string s = opt.extra_dm_internalprop_names[iextra]+to_string(opt.extra_dm_internalprop_index[iextra]);
            if (unique_extra_dm_internalprop_names.count(s)==0) {
                unique_extra_dm_internalprop_names.insert(s);
                opt.extra_dm_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_names_aperture.size();iextra++) {
            string s = opt.extra_dm_internalprop_names_aperture[iextra]+to_string(opt.extra_dm_internalprop_index_aperture[iextra]);
            if (unique_extra_dm_internalprop_names.count(s)==0) {
                unique_extra_dm_internalprop_names.insert(s);
                opt.extra_dm_internalprop_unique_input_indexlist.push_back(iextra);
            }
        }
        opt.extra_dm_internalprop_unique_input_names.assign(unique_extra_dm_internalprop_names.begin(), unique_extra_dm_internalprop_names.end());
#endif
}

///reads an hdf5 formatted file.
void ReadHDF(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons)
{
    //structure stores the names of the groups in the hdf input
    char buf[2000];
    HDF_Group_Names hdf_gnames (opt.ihdfnameconvention);
    //structures store names in groups
    vector<HDF_Header> hdf_header_info;
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
    vector<hid_t> Fhdf;
    vector<hid_t> partsgroup;
    vector<hid_t> headerattribs;
    vector<hid_t> headerdataspace;
    vector<hid_t> partsdataset;
    vector<hid_t> partsdataspace;
    vector<hid_t> partsdataset_extra;
    vector<hid_t> partsdataspace_extra;
    hid_t chunkspace;
    hid_t plist_id = H5P_DEFAULT;
    unsigned long long chunksize=opt.inputbufsize;
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
    Int_t i,j,k;
    unsigned long long n,nchunk,count,bcount,count2,bcount2;
    Int_t itemp;

    //store cosmology
    double z,aadjust,Hubble,Hubbleflow;

    double mscale,lscale,lvscale;
    double MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    int ifirstfile=0,*ireadfile,ireaderror=0;
    int *ireadtask,*readtaskID;
    Int_t ninputoffset;

    //for extra fields related to chemistry, feedback etc
    int numextrafields = 0;
    vector<int> numextrafieldsvec(NHDFTYPE);
    string extrafield, extrafield2;
    int iextraoffset;
    double *extrafieldbuff = NULL;

    SetUniqueInputNames(opt);

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
    //for parallel hdf5 read
    MPI_Comm mpi_comm_parallel_read;
    int ThisParallelReadTask, NProcsParallelReadTask;
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
    vector<hid_t> partsdatasetall;
    vector<hid_t> partsdataspaceall;
    vector<hid_t> partsdatasetall_extra;
    vector<hid_t> partsdataspaceall_extra;

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
        Nreadbuf=new Int_t[opt.nsnapread];
        for (int j=0;j<opt.nsnapread;j++) Nreadbuf[j]=0;
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
#ifdef USEPARALLELHDF
        if (opt.nsnapread > opt.num_files) {
            int ntaskread = ceil(opt.nsnapread/opt.num_files);
            int ifile = floor(ireadtask[ThisTask]/ntaskread);
	    int ThisReadTask, NProcsReadTask;
            MPI_Comm_rank(mpi_comm_read, &ThisReadTask);
            MPI_Comm_size(mpi_comm_read, &NProcsReadTask);
            MPI_Comm_split(mpi_comm_read, ifile, ThisReadTask, &mpi_comm_parallel_read);
            MPI_Comm_rank(mpi_comm_parallel_read, &ThisParallelReadTask);
            MPI_Comm_size(mpi_comm_parallel_read, &NProcsParallelReadTask);
            plist_id = H5Pcreate(H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(plist_id, mpi_comm_parallel_read, MPI_INFO_NULL);
        }
#endif
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
    hdf_header_info.resize(opt.num_files);
    for (i=0; i<opt.num_files; i++) hdf_header_info[i] = HDF_Header(opt.ihdfnameconvention);
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

    //handle any extra fields that should be loaded related to chemistry
    numextrafields = 0;
    for (auto &nf:numextrafieldsvec) nf=0;
    #if defined(GASON)
    numextrafieldsvec[HDFGASTYPE] = opt.gas_internalprop_names.size() + opt.gas_chem_names.size() + opt.gas_chemproduction_names.size();
    numextrafields += numextrafieldsvec[HDFGASTYPE];
    #endif
    #if defined(STARON)
    numextrafieldsvec[HDFSTARTYPE] = opt.star_internalprop_names.size() + opt.star_chem_names.size() + opt.star_chemproduction_names.size();
    numextrafields += numextrafieldsvec[HDFSTARTYPE];
    #endif
    #if defined(BHON)
    numextrafieldsvec[HDFBHTYPE] = opt.bh_internalprop_names.size() + opt.bh_chem_names.size() + opt.bh_chemproduction_names.size();
    numextrafields += numextrafieldsvec[HDFBHTYPE];
    #endif
    #if defined(EXTRADMON)
    numextrafieldsvec[HDFDMTYPE] = opt.extra_dm_internalprop_names.size();
    numextrafields += numextrafieldsvec[HDFDMTYPE];
    #endif
    if (numextrafields>0) {
        partsdataset_extra.resize(opt.num_files*numextrafields);
        partsdataspace_extra.resize(opt.num_files*numextrafields);
        for (auto &x:partsdataset_extra) x=-1;
        for (auto &x:partsdataspace_extra) x=-1;
    }

#ifdef USEMPI
    partsdatasetall.resize(opt.num_files*NHDFTYPE*NHDFDATABLOCK);
    partsdataspaceall.resize(opt.num_files*NHDFTYPE*NHDFDATABLOCK);
    for (auto &x:partsdatasetall) x=-1;
    for (auto &x:partsdataspaceall) x=-1;
    if (numextrafields>0) {
        partsdatasetall_extra.resize(opt.num_files*numextrafields);
        partsdataspaceall_extra.resize(opt.num_files*numextrafields);
        for (auto &x:partsdatasetall_extra) x=-1;
        for (auto &x:partsdataspaceall_extra) x=-1;
        extrafieldbuff = new double[numextrafields*chunksize];
    }
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
            Fhdf[i] = H5Fopen(buf, H5F_ACC_RDONLY, plist_id);
#ifdef USEPARALLELHDF
            H5Pclose(plist_id);
#endif
            if (ThisTask==0 && i==0) {
                cout<<buf<<endl;
                cout<<"HDF file contains the following group structures "<<endl;
                //H5Literate(Fhdf[i].getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
                H5Literate(Fhdf[i], H5_INDEX_NAME, H5_ITER_INC, NULL, file_info, NULL);
                cout<<" Expecting "<<endl;
                for (j=0;j<NHDFTYPE+1;j++) cout<<hdf_gnames.names[j]<<endl;
            }

            /* Read the BoxSize */
            if (opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES)  {
                /* SWIFT can have non-cubic boxes; but for cosmological runs they will always be cubes.
                * This makes the BoxSize a vector attribute, with it containing three values, but they
                * will always be the same. */
                hdf_header_info[i].BoxSize = read_attribute_v<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IBoxSize])[0];
            } else {
                hdf_header_info[i].BoxSize = read_attribute<double>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].IBoxSize]);
            }

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
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES) {
      mscale=opt.massinputconversion;lscale=opt.lengthinputconversion*aadjust;lvscale=opt.lengthinputconversion*opt.a;
    }
    else {
      mscale=opt.massinputconversion/opt.h;lscale=opt.lengthinputconversion/opt.h*aadjust;lvscale=opt.lengthinputconversion/opt.h*opt.a;
    }

    //ignore hubble flow
    Hubbleflow=0.;
    Ntotal=0;
    for (j=0;j<NHDFTYPE;j++)
    {
      opt.numpart[j]=hdf_header_info[ifirstfile].npartTotal[j];
      Ntotal+=hdf_header_info[ifirstfile].npartTotal[j];
    }
    for (j=0;j<NHDFTYPE;j++)
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
              /* If we have baryon search on, but ask to only search the dark matter, we'll segfault here.
               * Better to gracefully exit here with some helpful information. */

              cout << "\n";
              cout << "You have ran with the particle search type (=2) as Dark Matter but have \n";
              cout << "left the baryon search type as something nonzero. You should  set the \n";
              cout << "Baryon_searchflag to 0 in your parameter file.\n";

              cerr << "Incompatible choice of parameter values. See stdout for more information.\n";

#ifdef USE_MPI
              MPI_Abort(MPI_COMM_WORLD, 1);
#else
              exit(1);
#endif
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
#ifdef HIGHRES
                    else if (k==HDFDM1TYPE) Part[count].SetType(DARKTYPE);
                    else if (k==HDFDM2TYPE) Part[count].SetType(DARKTYPE);
#endif
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
#ifdef NOMASS
                  if (k==HDFDMTYPE) opt.MassValue = doublebuff[0];
#endif
                  for (int nn=0;nn<nchunk;nn++) Part[count++].SetMass(doublebuff[nn]);
                }
              }
              else {
#ifdef NOMASS
		if (k==HDFDMTYPE) opt.MassValue = hdf_header_info[i].mass[k];
#endif
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
                    for (int nn=0;nn<nchunk;nn++) Part[count++].SetSFR(doublebuff[nn] > 0. ? doublebuff[nn] : 0.);
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
                      for (int nn=0;nn<nchunk;nn++) Pbaryons[bcount++].SetZmet(doublebuff[nn]*zmetconversion);
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
            //now load extra fields if necessary.
            if (numextrafields>0)
            {
#if defined(GASON)
                if (opt.gas_internalprop_names.size() + opt.gas_chem_names.size() +
                    opt.gas_chemproduction_names.size()>0) {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k != HDFGASTYPE) {
                            count+=hdf_header_info[i].npart[k];
                            continue;
                        }
                        for(n=0;n<hdf_header_info[i].npart[k];n++) {
                            Part[count++].InitHydroProperties();
                        }
                    }
                }
                if (opt.gas_internalprop_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFGASTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.gas_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.gas_internalprop_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.gas_internalprop_names[iextra];
                                extrafield2 = extrafield + to_string(opt.gas_internalprop_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.gas_internalprop_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetHydroProperties().SetInternalProperties(extrafield2,doublebuff[nn]);
                                }
                                //close data spaces
                                for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                                for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                            }
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }

                if (opt.gas_chem_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFGASTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++)
                            for (auto &iextra:opt.gas_chem_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.gas_chem_names[iextra];
                                extrafield2 = extrafield + to_string(opt.gas_chem_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.gas_chem_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetHydroProperties().SetChemistry(extrafield2,doublebuff[nn]);
                                }
                                //close data spaces
                                for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                                for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                            }
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
                if (opt.gas_chemproduction_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFGASTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.gas_chemproduction_names.size();iextra++)
                            for (auto &iextra:opt.gas_chemproduction_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.gas_chemproduction_names[iextra];
                                extrafield2 = extrafield + to_string(opt.gas_chemproduction_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.gas_chemproduction_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetHydroProperties().SetChemistryProduction(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
#endif
#if defined(STARON)
                if (opt.star_internalprop_names.size() + opt.star_chem_names.size() +
                    opt.star_chemproduction_names.size()>0) {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k != HDFSTARTYPE) {
                            count+=hdf_header_info[i].npart[k];
                            continue;
                        }
                        for(n=0;n<hdf_header_info[i].npart[k];n++) {
                            Part[count++].InitStarProperties();
                        }
                    }
                }
                if (opt.star_internalprop_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFSTARTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.star_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.star_internalprop_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.star_internalprop_names[iextra];
                                extrafield2 = extrafield + to_string(opt.star_internalprop_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.star_internalprop_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetStarProperties().SetInternalProperties(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
                if (opt.star_chem_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFSTARTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.star_chem_names.size();iextra++)
                            for (auto &iextra:opt.star_chem_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.star_chem_names[iextra];
                                extrafield2 = extrafield + to_string(opt.star_chem_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.star_chem_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetStarProperties().SetChemistry(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
                if (opt.star_chemproduction_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFSTARTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.star_chemproduction_names.size();iextra++)
                            for (auto &iextra:opt.star_chemproduction_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.star_chemproduction_names[iextra];
                                extrafield2 = extrafield + to_string(opt.star_chemproduction_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.star_chemproduction_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetStarProperties().SetChemistryProduction(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
#endif
#if defined(BHON)
                if (opt.bh_internalprop_names.size() + opt.bh_chem_names.size() +
                    opt.bh_chemproduction_names.size()>0) {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k != HDFBHTYPE) {
                            count+=hdf_header_info[i].npart[k];
                            continue;
                        }
                        for(n=0;n<hdf_header_info[i].npart[k];n++) {
                            Part[count++].InitBHProperties();
                        }
                    }
                }
                if (opt.bh_internalprop_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFBHTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.bh_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.bh_internalprop_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.bh_internalprop_names[iextra];
                                extrafield2 = extrafield + to_string(opt.bh_internalprop_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.bh_internalprop_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetBHProperties().SetInternalProperties(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
                if (opt.bh_chem_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFBHTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.bh_chem_names.size();iextra++)
                            for (auto &iextra:opt.bh_chem_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.bh_chem_names[iextra];
                                extrafield2 = extrafield + to_string(opt.bh_chem_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.bh_chem_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetBHProperties().SetChemistry(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
                if (opt.bh_chemproduction_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFBHTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.bh_chemproduction_names.size();iextra++)
                            for (auto &iextra:opt.bh_chemproduction_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.bh_chemproduction_names[iextra];
                                extrafield2 = extrafield + to_string(opt.bh_chemproduction_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.bh_chemproduction_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetBHProperties().SetChemistryProduction(extrafield2,doublebuff[nn]);
                                }
                            }
                            //close data spaces
                            for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                            for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
#endif
#if defined(EXTRADMON)
                if (opt.extra_dm_internalprop_names.size()>0) {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k != HDFDMTYPE) {
                            count+=hdf_header_info[i].npart[k];
                            continue;
                        }
                        for(n=0;n<hdf_header_info[i].npart[k];n++) {
                            Part[count++].InitExtraDMProperties();
                        }
                    }
                }
                if (opt.extra_dm_internalprop_names.size()>0)
                {
                    count=count2;
                    bcount=bcount2;
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k == HDFDMTYPE)
                        {
                            // for (auto iextra=0;iextra<opt.extra_dm_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.extra_dm_internalprop_unique_input_indexlist)
                            {
                                unsigned long long count3=count;
                                extrafield = opt.extra_dm_internalprop_names[iextra];
                                extrafield2 = extrafield + to_string(opt.extra_dm_internalprop_index[iextra]);
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdataset_extra[i*numextrafields+iextra] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspace_extra[i*numextrafields+iextra] = HDF5OpenDataSpace(partsdataset_extra[i*numextrafields+iextra]);
                                //data loaded into memory in chunks
                                if (hdf_header_info[i].npart[k]<chunksize)nchunk = hdf_header_info[i].npart[k];
                                else nchunk=chunksize;
                                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                                {
                                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                                    //setup hyperslab so that it is loaded into the buffer
                                    HDF5ReadHyperSlabReal(doublebuff,partsdataset_extra[i*numextrafields+iextra],
                                        partsdataspace_extra[i*numextrafields+iextra], 1, 1, nchunk, n,
                                        plist_id, 1, opt.extra_dm_internalprop_index[iextra]);
                                    for (int nn=0;nn<nchunk;nn++) Part[count3++].GetExtraDMProperties().SetExtraProperties(extrafield2,doublebuff[nn]);
                                }
                                //close data spaces
                                for (auto &hidval:partsdataspace_extra) HDF5CloseDataSpace(hidval);
                                for (auto &hidval:partsdataset_extra) HDF5CloseDataSet(hidval);
                            }
                        }
                        else {
                            count+=hdf_header_info[i].npart[k];
                        }
                    }
                }
#endif
            }

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
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES) vscale = opt.velocityinputconversion;
    else vscale = opt.velocityinputconversion*sqrt(opt.a);
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES) opt.internalenergyinputconversion = opt.a*opt.a*opt.velocityinputconversion*opt.velocityinputconversion;
    else opt.internalenergyinputconversion = opt.velocityinputconversion*opt.velocityinputconversion;

    //finally adjust to appropriate units
    for (i=0;i<nbodies;i++)
    {
#ifdef HIGHRES
      if (Part[i].GetType()==DARKTYPE && Part[i].GetMass()<MP_DM) MP_DM=Part[i].GetMass();
#endif
      Part[i].SetMass(Part[i].GetMass()*mscale);
      for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*vscale+Hubbleflow*Part[i].GetPosition(j));
      for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
      for (i=0;i<nbaryons;i++)
      {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*vscale+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
      }
    }
#ifdef NOMASS
    opt.MassValue *= mscale;
#endif

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

                if (numextrafields>0)
                {
                    iextraoffset = 0;
#if defined(GASON)
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k!=HDFGASTYPE) continue;
                        if (opt.gas_internalprop_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.gas_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.gas_internalprop_unique_input_indexlist)
                            {
                                extrafield = opt.gas_internalprop_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.gas_internalprop_names.size();
                        if (opt.gas_chem_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++)
                            for (auto &iextra:opt.gas_chem_unique_input_indexlist)
                            {
                                extrafield = opt.gas_chem_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.gas_chem_names.size();
                        if (opt.gas_chemproduction_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.gas_chemproduction_names.size();iextra++)
                            for (auto &iextra:opt.gas_chemproduction_unique_input_indexlist)
                            {
                                extrafield = opt.gas_chemproduction_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.gas_chemproduction_names.size();
                    }
#endif
#if defined(STARON)
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k!=HDFSTARTYPE) continue;
                        if (opt.star_internalprop_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.star_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.star_internalprop_unique_input_indexlist)
                            {
                                extrafield = opt.star_internalprop_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.star_internalprop_names.size();
                        if (opt.star_chem_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.star_chem_names.size();iextra++)
                            for (auto &iextra:opt.star_chem_unique_input_indexlist)
                            {
                                extrafield = opt.star_chem_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.star_chem_names.size();
                        if (opt.star_chemproduction_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.star_chemproduction_names.size();iextra++)
                            for (auto &iextra:opt.star_chemproduction_unique_input_indexlist)
                            {
                                extrafield = opt.star_chemproduction_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.star_chemproduction_names.size();
                    }
#endif
#if defined(BHON)
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k!=HDFBHTYPE) continue;
                        if (opt.bh_internalprop_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.bh_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.bh_internalprop_unique_input_indexlist)
                            {
                                extrafield = opt.bh_internalprop_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.bh_internalprop_names.size();
                        if (opt.bh_chem_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.bh_chem_names.size();iextra++)
                            for (auto &iextra:opt.bh_chem_unique_input_indexlist)
                            {
                                extrafield = opt.bh_chem_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.bh_chem_names.size();
                        if (opt.bh_chemproduction_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.bh_chemproduction_names.size();iextra++)
                            for (auto &iextra:opt.bh_chemproduction_unique_input_indexlist)
                            {
                                extrafield = opt.bh_chemproduction_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.bh_chemproduction_names.size();
                    }
#endif
#if defined(EXTRADMON)
                    for (j=0;j<nusetypes;j++)
                    {
                        k=usetypes[j];
                        if (k!=HDFDMTYPE) continue;
                        if (opt.extra_dm_internalprop_names.size()>0)
                        {
                            // for (auto iextra=0;iextra<opt.extra_dm_internalprop_names.size();iextra++)
                            for (auto &iextra:opt.extra_dm_internalprop_unique_input_indexlist)
                            {
                                extrafield = opt.extra_dm_internalprop_names[iextra];
                                if (ThisTask==0 && opt.iverbose>1) cout<<"Opening group "<<hdf_gnames.part_names[k]<<": Data set "<<extrafield<<endl;
                                partsdatasetall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],extrafield);
                                partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset] = HDF5OpenDataSpace(partsdatasetall_extra[i*numextrafields+iextra+iextraoffset]);
                            }
                        }
                        iextraoffset += opt.extra_dm_internalprop_names.size();
                    }
#endif
                }

#ifdef USEPARALLELHDF
                if (opt.num_files<opt.nsnapread) {
                    plist_id = H5Pcreate(H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
                }
#endif
                for (j=0;j<nusetypes;j++)
                {
                    k=usetypes[j];
                    unsigned long long nstart = 0, nend = hdf_header_info[i].npart[k];
                    unsigned long long nlocalsize;
#ifdef USEPARALLELHDF
                    if (opt.num_files<opt.nsnapread) {
                        nlocalsize = nend / NProcsParallelReadTask;
                        nstart = nlocalsize*ThisParallelReadTask;
                        if (ThisParallelReadTask < NProcsParallelReadTask -1)
                            nend = nlocalsize + nstart;
                    }
#endif
                    if (nend-nstart<chunksize)nchunk=nend-nstart;
                    else nchunk=chunksize;
                    ninputoffset = 0;
                    for(n=nstart;n<nend;n+=nchunk)
                    {
                        if (nend - n < chunksize && nend - n > 0) nchunk=nend-n;
                        //setup hyperslab so that it is loaded into the buffer
                        //load positions
                        itemp=0;
                        //set hyperslab
                        HDF5ReadHyperSlabReal(doublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n, plist_id);
                        //velocities
                        itemp++;
                        HDF5ReadHyperSlabReal(veldoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n, plist_id);
                        //ids
                        itemp++;
                        HDF5ReadHyperSlabInteger(longbuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);

                        //masses
                        itemp++;
                        if (hdf_header_info[i].mass[k]==0) {
                            HDF5ReadHyperSlabReal(massdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                        }
#ifdef GASON
                        //self-energy
                        itemp++;
                        if (k == HDFGASTYPE) {
                            HDF5ReadHyperSlabReal(udoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                        }
#ifdef STARON
                        //star formation rate
                        itemp++;
                        if (k == HDFGASTYPE) {
                            HDF5ReadHyperSlabReal(SFRdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                        }

                        //metallicity
                        itemp++;
                        if (k == HDFGASTYPE || k == HDFSTARTYPE) {
                            HDF5ReadHyperSlabReal(Zdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                        }

                        //stellar age
                        itemp++;
                        if (k == HDFSTARTYPE) {
                            HDF5ReadHyperSlabReal(Tagedoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                        }
#endif
#endif
                        //load extra fields
                        if (numextrafields>0)
                        {
                            iextraoffset = 0;
#if defined(GASON)
                            if (opt.gas_internalprop_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.gas_internalprop_names.size();iextra++)
                                for (auto &iextra:opt.gas_internalprop_unique_input_indexlist)
                                {
                                    if (k == HDFGASTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.gas_internalprop_index[iextra]);
                                }
                            }
                            iextraoffset += opt.gas_internalprop_names.size();
                            if (opt.gas_chem_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++)
                                for (auto &iextra:opt.gas_chem_unique_input_indexlist)
                                {
                                    if (k == HDFGASTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.gas_chem_index[iextra]);
                                }
                            }
                            iextraoffset += opt.gas_chem_names.size();
                            if (opt.gas_chemproduction_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.gas_chemproduction_names.size();iextra++)
                                for (auto &iextra:opt.gas_chemproduction_unique_input_indexlist)
                                {
                                    if (k == HDFGASTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.gas_chemproduction_index[iextra]);
                                }
                            }
                            iextraoffset += opt.gas_chemproduction_names.size();
#endif
#if defined(STARON)
                            if (opt.star_internalprop_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.star_internalprop_names.size();iextra++)
                                for (auto &iextra:opt.star_internalprop_unique_input_indexlist)
                                {
                                    if (k == HDFSTARTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.star_internalprop_index[iextra]);
                                }
                            }
                            iextraoffset += opt.star_internalprop_names.size();
                            if (opt.star_chem_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.star_chem_names.size();iextra++)
                                for (auto &iextra:opt.star_chem_unique_input_indexlist)
                                {
                                    if (k == HDFSTARTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.star_chem_index[iextra]);
                                }
                            }
                            iextraoffset += opt.star_chem_names.size();
                            if (opt.star_chemproduction_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.star_chemproduction_names.size();iextra++)
                                for (auto &iextra:opt.star_chemproduction_unique_input_indexlist)
                                {
                                    if (k == HDFSTARTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.star_chemproduction_index[iextra]);
                                }
                            }
                            iextraoffset += opt.star_chemproduction_names.size();
#endif
#if defined(BHON)
                            if (opt.bh_internalprop_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.bh_internalprop_names.size();iextra++)
                                for (auto &iextra:opt.bh_internalprop_unique_input_indexlist)
                                {
                                    if (k == HDFBHTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.bh_internalprop_index[iextra]);
                                }
                            }
                            iextraoffset += opt.bh_internalprop_names.size();
                            if (opt.bh_chem_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.bh_chem_names.size();iextra++)
                                for (auto &iextra:opt.bh_chem_unique_input_indexlist)
                                {
                                    if (k == HDFBHTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.bh_chem_index[iextra]);
                                }
                            }
                            iextraoffset += opt.bh_chem_names.size();
                            if (opt.bh_chemproduction_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.bh_chemproduction_names.size();iextra++)
                                for (auto &iextra:opt.bh_chemproduction_unique_input_indexlist)
                                {
                                    if (k == HDFBHTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.bh_chemproduction_index[iextra]);
                                }
                            }
                            iextraoffset += opt.bh_chemproduction_names.size();
#endif
#if defined(EXTRADMON)
                            if (opt.extra_dm_internalprop_names.size()>0)
                            {
                                // for (auto iextra=0;iextra<opt.extra_dm_internalprop_names.size();iextra++)
                                for (auto &iextra:opt.extra_dm_internalprop_unique_input_indexlist)
                                {
                                    if (k == HDFDMTYPE)
                                        HDF5ReadHyperSlabReal(&extrafieldbuff[(iextraoffset+iextra)*chunksize],partsdatasetall_extra[i*numextrafields+iextra+iextraoffset], partsdataspaceall_extra[i*numextrafields+iextra+iextraoffset], 1, 1, nchunk, n, plist_id, 1, opt.extra_dm_internalprop_index[iextra]);
                                }
                            }
                            iextraoffset += opt.extra_dm_internalprop_names.size();
#endif
                        }
                    for (unsigned long long nn=0;nn<nchunk;nn++) {
                        ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        ibufindex=ibuf*BufSize+Nbuf[ibuf];
                        //reset hydro quantities of buffer
#ifdef GASON
                        Pbuf[ibufindex].SetU(0);
                        Pbuf[ibufindex].SetHydroProperties();
#ifdef STARON
                        Pbuf[ibufindex].SetSFR(0);
                        Pbuf[ibufindex].SetZmet(0);
#endif
#endif
#ifdef STARON
                        Pbuf[ibufindex].SetZmet(0);
                        Pbuf[ibufindex].SetTage(0);
                        Pbuf[ibufindex].SetStarProperties();
#endif
#ifdef BHON
                        Pbuf[ibufindex].SetBHProperties();
#endif
#ifdef EXTRADMON
                        Pbuf[ibufindex].SetExtraDMProperties();
#endif

                        Pbuf[ibufindex].SetPosition(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Pbuf[ibufindex].SetVelocity(veldoublebuff[nn*3],veldoublebuff[nn*3+1],veldoublebuff[nn*3+2]);
#ifdef NOMASS
                        if (k==HDFDMTYPE) {
                            if (hdf_header_info[i].mass[k]==0) opt.MassValue = massdoublebuff[0];
                            else opt.MassValue = hdf_header_info[i].mass[k];
                        }
#endif
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
                        if (k==HDFGASTYPE && MP_B<Pbuf[ibufindex].GetMass()) MP_B=Pbuf[ibufindex].GetMass();
#endif
#ifdef GASON
                      if (k==HDFGASTYPE) {
                        Pbuf[ibufindex].SetU(udoublebuff[nn]);
#ifdef STARON
                        Pbuf[ibufindex].SetSFR(SFRdoublebuff[nn] > 0. ? SFRdoublebuff[nn]: 0.);
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

                    if (numextrafields>0) {
                        iextraoffset = 0;
#ifdef GASON
                        if (k==HDFGASTYPE && numextrafieldsvec[HDFGASTYPE]) {
                            if (!Pbuf[ibufindex].HasHydroProperties()) Pbuf[ibufindex].InitHydroProperties();
                            if (opt.gas_internalprop_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.gas_internalprop_names.size();iextra++)
                                {
                                    extrafield = opt.gas_internalprop_names[iextra] +
                                        to_string(opt.gas_internalprop_index[iextra]);
                                    Pbuf[ibufindex].GetHydroProperties().SetInternalProperties(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.gas_internalprop_names.size();
                            if (opt.gas_chem_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++)
                                {
                                    extrafield = opt.gas_chem_names[iextra] +
                                        to_string(opt.gas_chem_index[iextra]);
                                    Pbuf[ibufindex].GetHydroProperties().SetChemistry(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.gas_chem_names.size();
                            if (opt.gas_chemproduction_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.gas_chemproduction_names.size();iextra++)
                                {
                                    extrafield = opt.gas_chemproduction_names[iextra] +
                                        to_string(opt.gas_chemproduction_index[iextra]);
                                    Pbuf[ibufindex].GetHydroProperties().SetChemistryProduction(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.gas_chemproduction_names.size();
                        }
#endif
#ifdef STARON
                        if (k==HDFSTARTYPE && numextrafieldsvec[HDFSTARTYPE]) {
                            if (!Pbuf[ibufindex].HasStarProperties()) Pbuf[ibufindex].InitStarProperties();
                            if (opt.star_internalprop_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.star_internalprop_names.size();iextra++)
                                {
                                    extrafield = opt.star_internalprop_names[iextra] +
                                        to_string(opt.star_internalprop_index[iextra]);
                                    Pbuf[ibufindex].GetStarProperties().SetInternalProperties(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.star_internalprop_names.size();
                            if (opt.star_chem_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.star_chem_names.size();iextra++)
                                {
                                    extrafield = opt.star_chem_names[iextra] +
                                        to_string(opt.star_chem_index[iextra]);
                                    Pbuf[ibufindex].GetStarProperties().SetChemistry(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.star_chem_names.size();
                            if (opt.star_chemproduction_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.star_chemproduction_names.size();iextra++)
                                {
                                    extrafield = opt.star_chemproduction_names[iextra] +
                                        to_string(opt.star_chemproduction_index[iextra]);
                                    Pbuf[ibufindex].GetStarProperties().SetChemistryProduction(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.star_chemproduction_names.size();
                        }
#endif
#ifdef BHON
                        if (k==HDFBHTYPE && numextrafieldsvec[HDFBHTYPE]) {
                            if (!Pbuf[ibufindex].HasBHProperties()) Pbuf[ibufindex].InitBHProperties();
                            if (opt.bh_internalprop_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.bh_internalprop_names.size();iextra++)
                                {
                                    extrafield = opt.bh_internalprop_names[iextra] +
                                        to_string(opt.bh_internalprop_index[iextra]);
                                    Pbuf[ibufindex].GetBHProperties().SetInternalProperties(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.bh_internalprop_names.size();
                            if (opt.bh_chem_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.bh_chem_names.size();iextra++)
                                {
                                    extrafield = opt.bh_chem_names[iextra] +
                                        to_string(opt.bh_chem_index[iextra]);
                                    Pbuf[ibufindex].GetBHProperties().SetChemistry(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.bh_chem_names.size();
                            if (opt.bh_chemproduction_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.bh_chemproduction_names.size();iextra++)
                                {
                                    extrafield = opt.bh_chemproduction_names[iextra] +
                                        to_string(opt.bh_chemproduction_index[iextra]);
                                    Pbuf[ibufindex].GetBHProperties().SetChemistryProduction(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.bh_chemproduction_names.size();
                        }
#endif
#ifdef EXTRADMON
                        if (k==HDFDMTYPE && numextrafieldsvec[HDFDMTYPE]) {
                            if (!Pbuf[ibufindex].HasExtraDMProperties()) Pbuf[ibufindex].InitExtraDMProperties();
                            if (opt.extra_dm_internalprop_names.size()>0)
                            {
                                for (auto iextra=0;iextra<opt.extra_dm_internalprop_names.size();iextra++)
                                {
                                    extrafield = opt.extra_dm_internalprop_names[iextra] +
                                        to_string(opt.extra_dm_internalprop_index[iextra]);
                                    Pbuf[ibufindex].GetExtraDMProperties().SetExtraProperties(extrafield, extrafieldbuff[(iextraoffset+iextra)*chunksize+nn]);
                                }
                            }
                            iextraoffset += opt.extra_dm_internalprop_names.size();
                        }
#endif
                    }

#ifdef EXTRAINPUTINFO
                        if (opt.iextendedoutput)
                        {
                            Pbuf[ibufindex].SetInputFileID(i);
                            Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                        }
#endif
                      Nbuf[ibuf]++;
                      MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
                    }
                    ninputoffset += nchunk;
                  }
                }
                if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                  for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    unsigned long long nstart = 0, nend = hdf_header_info[i].npart[k];
                    unsigned long long nlocalsize;
#ifdef USEPARALLELHDF
                    if (opt.num_files<opt.nsnapread) {
                        nlocalsize = nend / NProcsParallelReadTask;
                        nstart = nlocalsize*ThisParallelReadTask;
                        if (ThisParallelReadTask < NProcsParallelReadTask -1)
                            nend = nlocalsize + nstart;
                    }
#endif
                    if (nend-nstart<chunksize)nchunk=nend-nstart;
                    else nchunk=chunksize;
                    ninputoffset = 0;
                    for(n=nstart;n<nend;n+=nchunk)
                    {
                      if (nend - n < chunksize && nend - n > 0) nchunk=nend-n;
                      //setup hyperslab so that it is loaded into the buffer
                      //load positions
                      itemp=0;
                      //set hyperslab
                      HDF5ReadHyperSlabReal(doublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n, plist_id);
                      //velocities
                      itemp++;
                      HDF5ReadHyperSlabReal(veldoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n, plist_id);
                      //ids
                      itemp++;
                      HDF5ReadHyperSlabInteger(longbuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 3, nchunk, n, plist_id);
                      //masses
                      itemp++;
                      if (hdf_header_info[i].mass[k]==0) {
                          HDF5ReadHyperSlabReal(massdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                      }
#ifdef GASON
                      //self-energy
                      itemp++;
                      if (k == HDFGASTYPE) {
                          HDF5ReadHyperSlabReal(udoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                      }
#ifdef STARON
                      //star formation rate
                      itemp++;
                      if (k == HDFGASTYPE) {
                          HDF5ReadHyperSlabReal(SFRdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                      }

                      //metallicity
                      itemp++;
                      if (k == HDFGASTYPE || k == HDFSTARTYPE) {
                          HDF5ReadHyperSlabReal(Zdoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
                      }

                      //stellar age
                      itemp++;
                      if (k == HDFSTARTYPE) {
                          HDF5ReadHyperSlabReal(Tagedoublebuff,partsdatasetall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], partsdataspaceall[i*NHDFTYPE*NHDFDATABLOCK+k*NHDFDATABLOCK+itemp], 1, 1, nchunk, n, plist_id);
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
                          Pbuf[ibufindex].SetSFR(SFRdoublebuff[nn] > 0. ? SFRdoublebuff[nn] : 0.);
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
                        MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
                      }
                      ninputoffset+=nchunk;
                    }//end of chunk
                  }//end of part type
                }//end of baryon if
                //close property
#ifdef USEPARALLELHDF
                H5Pclose(plist_id);
#endif
                //close data spaces
                for (auto &hidval:partsdataspaceall) HDF5CloseDataSpace(hidval);
                for (auto &hidval:partsdatasetall) HDF5CloseDataSet(hidval);
                for (auto &hidval:partsdataspaceall_extra) HDF5CloseDataSpace(hidval);
                for (auto &hidval:partsdatasetall_extra) HDF5CloseDataSet(hidval);
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
                MPISendHydroInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
                MPISendStarInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
                MPISendBHInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
                MPISendExtraDMInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
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
      if (opt.Omega_b==0 && MP_B>0){
        opt.Omega_b=MP_B/(MP_DM+MP_B)*opt.Omega_m;
        opt.Omega_cdm=opt.Omega_m-opt.Omega_b;
      }

      // SWIFT snapshots already include the 1/h factor factor,
      // so there is no need to include it.
      if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES)
      {
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
#ifdef HIGHRES
    if (opt.nsnapread>1) {
        MPI_Allreduce(MPI_IN_PLACE,&MP_DM, 1, MPI_DOUBLE, MPI_MIN,mpi_comm_read);
    }
#endif
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
    opt.MassValue *= mscale;
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
#ifdef USEPARALLELHDF
    if (opt.nsnapread > opt.num_files) MPI_Comm_free(&mpi_comm_parallel_read);
#endif
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
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES) vscale = opt.velocityinputconversion;
    else vscale = opt.velocityinputconversion*sqrt(opt.a);
    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES) opt.internalenergyinputconversion = opt.a*opt.a*opt.velocityinputconversion*opt.velocityinputconversion;
    else opt.internalenergyinputconversion = opt.velocityinputconversion*opt.velocityinputconversion;


    //finally adjust to appropriate units
    for (i=0;i<Nlocal;i++)
    {
      Part[i].SetMass(Part[i].GetMass()*mscale);
      for (int j=0;j<3;j++) Part[i].SetVelocity(j,Part[i].GetVelocity(j)*vscale+Hubbleflow*Part[i].GetPosition(j));
      for (int j=0;j<3;j++) Part[i].SetPosition(j,Part[i].GetPosition(j)*lscale);
    }
    if (Pbaryons!=NULL && opt.iBaryonSearch==1) {
      for (i=0;i<Nlocalbaryon[0];i++)
      {
        Pbaryons[i].SetMass(Pbaryons[i].GetMass()*mscale);
        for (int j=0;j<3;j++) Pbaryons[i].SetVelocity(j,Pbaryons[i].GetVelocity(j)*vscale+Hubbleflow*Pbaryons[i].GetPosition(j));
        for (int j=0;j<3;j++) Pbaryons[i].SetPosition(j,Pbaryons[i].GetPosition(j)*lscale);
    }
    }
#endif

    delete[] intbuff;
    delete[] longbuff;
    delete[] uintbuff;
    delete[] floatbuff;
    delete[] doublebuff;
    delete[] extrafieldbuff;
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

    UpdateExtraFieldNames(opt);
}

#endif
