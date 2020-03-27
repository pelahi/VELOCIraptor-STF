/*! \file hdfitems.h
 *  \brief this file contains definitions and routines for reading HDF5 files
 *
 *   NOTE: the routines are based on reading Illustris HDF outputs
 */

#ifndef HDFITEMS_H
#define HDFITEMS_H

#include "hdf5.h"


///\name ILLUSTRIS specific constants
//@{
///convert illustris metallicty to ratio to solar
#define ILLUSTRISZMET 1.0/0.0127
//@}

///number of particle types
#define NHDFTYPE 6
///\name define hdf particle types
//@{
#define HDFGASTYPE 0
#define HDFDMTYPE 1
#define HDFWINDTYPE 2
#define HDFTRACERTYPE 3
#define HDFSTARTYPE 4
#define HDFBHTYPE 5
#define HDFEXTRATYPE 6

#define HDFDM1TYPE 2
#define HDFDM2TYPE 3
//@}

///\name number of entries in various data groups
//@{
#define HDFHEADNINFO 12
#define HDFGASNINFO 20
#define HDFDMNINFO 7
#define HDFTRACERNINFO 3
#define HDFSTARNINFO 13
#define HDFBHNINFO 21
#define HDFMAXPINFO 40
//@}

///\name to access where extra baryonic properties are located in the HDF_Part_Info structure that code will use to calcualte object properties
//@{
#define HDFGASIMETAL 0

#define HDFSTARIMETAL 40
#define HDFSTARIAGE 41

#define HDFBHIMETAL 50
#define HDFBHIAGE 51
#define HDFBHIMDOT 52
//@}

///number of luminosity bands for stars
#define NLUMBANDS 8

///how many particle properties are read from file in one go
#define HDFCHUNKSIZE 100000

///how many data blocks are of  interest. That is what data blocks do I wish to load. pos,vel,mass,pid,U,Zmet,tage?
#define NHDFDATABLOCK 10
///here number shared by all particle types
#define NHDFDATABLOCKALL 4
//Maximum dimensionality of a datablock
///example at most one needs a dimensionality of 13 for the tracer particles in Illustris for fluid related info
#define HDFMAXPROPDIM 13

///\defgroup HDFNAMES labels for HDF naming conventions
//@{
#define HDFNUMNAMETYPES  9
#define HDFILLUSTISNAMES 0
#define HDFGADGETXNAMES  1
#define HDFEAGLENAMES    2
#define HDFGIZMONAMES    3
#define HDFSIMBANAMES    4
#define HDFMUFASANAMES   5
#define HDFSWIFTEAGLENAMES    6
#define HDFOLDSWIFTEAGLENAMES    8
#define HDFEAGLEVERSION2NAMES    7
//@}

///size of chunks in hdf files for Compression
#define HDFOUTPUTCHUNKSIZE 8192
#define HDFDEFLATE    6


#if H5_VERSION_GE(1,10,1)
#define HDF5_FILE_GROUP_COMMON_BASE H5::Group
#else
#define HDF5_FILE_GROUP_COMMON_BASE H5::CommonFG
#endif

// #ifdef USEPARALLELHDF
// // #if H5_VERSION_GE(1,10,2)
// // #define USEHDFCOMPRESSOIN
// // #endif
// // #else
// // #define USEHDFCOMPRESSOIN
// #endif

template <typename ReturnT, typename F, typename ... Ts>
ReturnT safe_hdf5(F function, Ts ... args)
{
       ReturnT status = function(std::forward<Ts>(args)...);
       if (status < 0) {
           cerr<<"Error in HDF routine "<<endl;//<<function.__PRETTY_FUNCTION__
           //throw std::runtime_error("Error in HDF routine.");
           #ifdef USEMPI
           MPI_Abort(MPI_COMM_WORLD,1);
           #else
           exit(1);
           #endif
       }
       return status;
}

// Overloaded function to return HDF5 type given a C type
static inline hid_t hdf5_type(float dummy)              {return H5T_NATIVE_FLOAT;}
static inline hid_t hdf5_type(double dummy)             {return H5T_NATIVE_DOUBLE;}
static inline hid_t hdf5_type(short dummy)              {return H5T_NATIVE_SHORT;}
static inline hid_t hdf5_type(int dummy)                {return H5T_NATIVE_INT;}
static inline hid_t hdf5_type(long dummy)               {return H5T_NATIVE_LONG;}
static inline hid_t hdf5_type(long long dummy)          {return H5T_NATIVE_LLONG;}
static inline hid_t hdf5_type(unsigned short dummy)     {return H5T_NATIVE_USHORT;}
static inline hid_t hdf5_type(unsigned int dummy)       {return H5T_NATIVE_UINT;}
static inline hid_t hdf5_type(unsigned long dummy)      {return H5T_NATIVE_ULONG;}
static inline hid_t hdf5_type(unsigned long long dummy) {return H5T_NATIVE_ULLONG;}
static inline hid_t hdf5_type(std::string dummy)        {return H5T_C_S1;}

//template <typename AttributeHolder>
//static inline H5::Attribute get_attribute(const AttributeHolder &l, const std::string attr_name)
static inline void get_attribute(vector<hid_t> &ids, const std::string attr_name)
{
    //can use H5Aexists as it is the C interface but how to access it?
    //auto exists = H5Aexists(l.getId(), attr_name.c_str());
    auto exists = H5Aexists(ids.back(), attr_name.c_str());
    if (exists == 0) {
        throw invalid_argument(std::string("attribute not found ") + attr_name);
    }
    else if (exists < 0) {
        throw std::runtime_error("Error on H5Aexists");
    }
    auto attr = H5Aopen(ids.back(), attr_name.c_str(), H5P_DEFAULT);
    ids.push_back(attr);
}


static inline void get_attribute(vector<hid_t> &ids, const std::vector<std::string> &parts)
{
    // This is the attribute name, so open it and store the id
    if (parts.size() == 1) {
        get_attribute(ids, parts[0]);
    }
    else {
        H5O_info_t object_info;
        hid_t newid;
        H5Oget_info_by_name(ids.back(), parts[0].c_str(), &object_info, H5P_DEFAULT);
        if (object_info.type == H5O_TYPE_GROUP) {
            newid = H5Gopen2(ids.back(),parts[0].c_str(),H5P_DEFAULT);
        }
        else if (object_info.type == H5O_TYPE_DATASET) {
            newid = H5Dopen2(ids.back(),parts[0].c_str(),H5P_DEFAULT);
        }
        ids.push_back(newid);
        //get the substring
        vector<string> subparts(parts.begin() + 1, parts.end());
        //call function again
        get_attribute(ids, subparts);
    }
    //throw invalid_argument("attribute name not found");
}

static inline vector<string> tokenize(const string &s, const string &delims)
{
    string::size_type lastPos = s.find_first_not_of(delims, 0);
    string::size_type pos     = s.find_first_of(delims, lastPos);

    vector<string> tokens;
    while (string::npos != pos || string::npos != lastPos) {
        tokens.push_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delims, pos);
        pos = s.find_first_of(delims, lastPos);
    }
    return tokens;
}

static inline void get_attribute(const hid_t &file_id, vector<hid_t> &ids, const string &name)
{
    std::vector<std::string> parts = tokenize(name, "/");
    ids.push_back(file_id);
    get_attribute(ids, parts);
}

static inline void close_hdf_ids(vector<hid_t> &ids)
{
    H5O_info_t object_info;
    for (auto &id:ids)
    {
        H5Oget_info(id, &object_info);
        if (object_info.type == H5O_TYPE_GROUP) {
            H5Gclose(id);
        }
        else if (object_info.type == H5O_TYPE_GROUP) {
            H5Dclose(id);
        }
    }
}

template<typename T> static inline void _do_read(const hid_t &attr, const hid_t &type, T &val)
{
    H5Aread(attr, type, &val);
}

template<> void _do_read<std::string>(const hid_t &attr, const hid_t &type, std::string &val)
{
    vector<char> buf;
    hid_t type_in_file = H5Aget_type(attr);
    hid_t type_in_memory = H5Tcopy(type); // copy memory type because we'll need to modify it
    size_t length = H5Tget_size(type_in_file); // get length of the string in the file
    buf.resize(length+1); // resize buffer in memory, allowing for null terminator
    H5Tset_size(type_in_memory, length+1); // tell HDF5 the length of the buffer in memory
    H5Tset_strpad(type_in_memory, H5T_STR_NULLTERM); // specify that we want a null terminated string
    H5Aread(attr, type_in_memory, buf.data());
    H5Tclose(type_in_memory);
    H5Tclose(type_in_file);
    val=string(buf.data());
}

template<typename T> static inline void _do_read_v(const hid_t &attr, const hid_t &type, vector<T> &val)
{
    hid_t space = H5Aget_space (attr);
    int npoints = H5Sget_simple_extent_npoints(space);
    val.resize(npoints);
    H5Aread(attr, type, val.data());
    H5Sclose(space);
}

template<typename T> const T read_attribute(const hid_t &file_id, const std::string &name) {
    std::string attr_name;
    T val;
    hid_t type;
    H5O_info_t object_info;
    vector <hid_t> ids;
    //traverse the file to get to the attribute, storing the ids of the
    //groups, data spaces, etc that have been opened.
    get_attribute(file_id, ids, name);
    //now reverse ids and load attribute
    reverse(ids.begin(),ids.end());
    //determine hdf5 type of the array in memory
    type = hdf5_type(T{});
    // read the data
    _do_read<T>(ids[0], type, val);
    H5Aclose(ids[0]);
    //remove file id from id list
    ids.pop_back();
    ids.erase(ids.begin());
    //now have hdf5 ids traversed to get to desired attribute so move along to close all
    //based on their object type
    close_hdf_ids(ids);
    return val;
}

//read vector attribute
template<typename T> const vector<T> read_attribute_v(const hid_t &file_id, const std::string &name) {
    std::string attr_name;
    vector<T> val;
    hid_t type;
    H5O_info_t object_info;
    vector <hid_t> ids;
    //traverse the file to get to the attribute, storing the ids of the
    //groups, data spaces, etc that have been opened.
    get_attribute(file_id, ids, name);
    //now reverse ids and load attribute
    reverse(ids.begin(),ids.end());
    //determine hdf5 type of the array in memory
    type = hdf5_type(T{});
    // read the data
    _do_read_v<T>(ids[0], type, val);
    H5Aclose(ids[0]);
    //remove file id from id list
    ids.pop_back();
    ids.erase(ids.begin());
    //now have hdf5 ids traversed to get to desired attribute so move along to close all
    //based on their object type
    close_hdf_ids(ids);
    return val;
}


template<typename T> const T read_attribute(const std::string &filename, const std::string &name) {
    safe_hdf5<herr_t>(H5Fopen, filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    T attr = read_attribute<T>(file_id, name);
    safe_hdf5<herr_t>(H5Fclose,file_id);
    return attr;
}

static inline hid_t HDF5OpenFile(string name, unsigned int flags){
    hid_t Fhdf;
    return H5Fopen(name.c_str(),flags, H5P_DEFAULT);
}

static inline hid_t HDF5OpenGroup(const hid_t &file, string name){
    return H5Gopen2(file,name.c_str(),H5P_DEFAULT);
}
static inline hid_t HDF5OpenDataSet(const hid_t &id, string name){
    hid_t idval = H5Dopen2(id,name.c_str(),H5P_DEFAULT);
    return idval;
}
static inline hid_t HDF5OpenDataSpace(const hid_t &id){
    hid_t idval=H5Dget_space(id);
    return idval;
}


static inline int whatisopen(hid_t fid) {
        ssize_t cnt;
        int howmany;
        int i;
        H5I_type_t ot;
        hid_t anobj;
        hid_t *objs;
        char name[1024];
        herr_t status;

        cnt = H5Fget_obj_count(fid, H5F_OBJ_ALL);

        if (cnt <= 0) return cnt;

        printf("%zd object(s) open\n", cnt);

        objs = new hid_t[cnt];

        howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, cnt, objs);

        printf("open objects:\n");

        for (i = 0; i < howmany; i++ ) {
             anobj = objs[i];
             ot = H5Iget_type(anobj);
             status = H5Iget_name(anobj, name, 1024);
             printf(" %d: type %d, name %s\n",i,ot,name);
        }
	delete[] objs;
        return howmany;
}


static inline void HDF5CloseFile(hid_t &id){
    if (id>=0) H5Fclose(id);
    id = -1;
}
static inline void HDF5CloseGroup(hid_t &id){
    if (id>=0) H5Gclose(id);
    id = -1;
}
static inline void HDF5CloseDataSet(hid_t &id){
    if (id>=0) H5Dclose(id);
    id = -1;
}
static inline void HDF5CloseDataSpace(hid_t &id){
    if (id>=0) H5Sclose(id);
    id = -1;
}

static inline void HDF5ReadHyperSlabReal(double *buffer,
    const hid_t &dataset, const hid_t &dataspace,
    const hsize_t datarank, const hsize_t dummy,
    unsigned long long nchunk, unsigned long long noffset,
    hid_t plist_id = H5P_DEFAULT,
    unsigned long long nchunk2 = 0, unsigned long long noffset2 = 0,
    vector<hsize_t> nchunkvec = vector<hsize_t>(),
    vector<hsize_t> noffsetvec = vector<hsize_t>()
)
{
    //setup hyperslab so that it is loaded into the buffer
    vector<hsize_t> start, count, stride, block, memdims, dims;
    hsize_t ndim, memsize = 1;
    hid_t memspace;

    //init for 1d array
    start.push_back(noffset);
    count.push_back(nchunk);
    stride.push_back(1);
    block.push_back(1);

    //get dimensions of data set
    ndim = H5Sget_simple_extent_ndims(dataspace);
    dims.resize(ndim);
    //get extent in each dimension
    H5Sget_simple_extent_dims(dataspace, dims.data(), NULL);
    //if offsets not explicilty provided
    //then assume offset is zero for any higher dimensions
    if (noffsetvec.size() == 0) {
        noffsetvec.resize(ndim,0);
    }
    //if chunksize vector not provided, assume extent
    if (nchunkvec.size() == 0) {
        nchunkvec = dims;
    }
    if (ndim==2) {
        if (nchunk2 == 0) nchunk2 = dims[1];
        start.push_back(noffset2);
        count.push_back(nchunk2);
        stride.push_back(1);
        block.push_back(1);
    }
    else if (ndim>2) {
        for (auto i=1;i<ndim;i++) {
            start.push_back(noffsetvec[i]);
            count.push_back(nchunkvec[i]);
            stride.push_back(1);
            block.push_back(1);
        }
    }
    for (auto x:count) memsize *= x;
    //if the return data rank is 1, then all data stored in 1d array
    if (datarank > 1) {
        //not implemented yet, still map everything to 1d array
        memdims.push_back(memsize);
    }
    else {
        memdims.push_back(memsize);
    }
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), stride.data(), count.data(), block.data());
    memspace = H5Screate_simple (1, memdims.data(), NULL);
    safe_hdf5<herr_t>(H5Dread, dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, buffer);

}

static inline void HDF5ReadHyperSlabInteger(long long *buffer,
    const hid_t &dataset, const hid_t &dataspace,
    const hsize_t datarank, const hsize_t dummy,
    unsigned long long nchunk, unsigned long long noffset,
    hid_t plist_id = H5P_DEFAULT,
    unsigned long long nchunk2 = 0, unsigned long long noffset2 = 0,
    vector<hsize_t> nchunkvec = vector<hsize_t>(),
    vector<hsize_t> noffsetvec = vector<hsize_t>()
)
{
    //setup hyperslab so that it is loaded into the buffer
    vector<hsize_t> start, count, stride, block, memdims, dims;
    hsize_t ndim, memsize = 1;
    hid_t memspace;

    //init for 1d array
    start.push_back(noffset);
    count.push_back(nchunk);
    stride.push_back(1);
    block.push_back(1);

    //get dimensions of data set
    ndim = H5Sget_simple_extent_ndims(dataspace);
    dims.resize(ndim);
    //get extent in each dimension
    H5Sget_simple_extent_dims(dataspace, dims.data(), NULL);
    //if offsets not explicilty provided
    //then assume offset is zero for any higher dimensions
    if (noffsetvec.size() == 0) {
        noffsetvec.resize(ndim,0);
    }
    //if chunksize vector not provided, assume extent
    if (nchunkvec.size() == 0) {
        nchunkvec = dims;
    }
    if (ndim==2) {
        if (nchunk2 == 0) nchunk2 = dims[1];
        start.push_back(nchunk2);
        count.push_back(noffset2);
        stride.push_back(1);
        block.push_back(1);
    }
    else if (ndim>2) {
        for (auto i=1;i<ndim;i++) {
            start.push_back(noffsetvec[i]);
            count.push_back(nchunkvec[i]);
            stride.push_back(1);
            block.push_back(1);
        }
    }
    for (auto x:count) memsize *= x;
    //if the return data rank is 1, then all data stored in 1d array
    if (datarank > 1) {
        //not implemented yet, still map everything to 1d array
        memdims.push_back(memsize);
    }
    else {
        memdims.push_back(memsize);
    }
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), stride.data(), count.data(), block.data());
    memspace = H5Screate_simple (1, memdims.data(), NULL);
    safe_hdf5<herr_t>(H5Dread, dataset, H5T_NATIVE_LONG, memspace, dataspace, plist_id, buffer);
}

///\name HDF class to manage writing information
class H5OutputFile
{
    protected:

    hid_t file_id;
#ifdef USEPARALLELHDF
    hid_t parallel_access_id;
#endif

    // Called if a HDF5 call fails (might need to MPI_Abort)
    void io_error(std::string message) {
    std::cerr << message << std::endl;
#ifdef USEMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    abort();
    }

    public:

    // Constructor
    H5OutputFile() {
        file_id = -1;
#ifdef USEPARALLELHDF
        parallel_access_id = -1;
#endif
    }

    // Create a new file
    void create(std::string filename, hid_t flag = H5F_ACC_TRUNC,
        int taskID = -1, bool iparallelopen = true)
    {
        if(file_id >= 0)io_error("Attempted to create file when already open!");
#ifdef USEPARALLELHDF
        MPI_Comm comm = mpi_comm_write;
        MPI_Info info = MPI_INFO_NULL;
        if (iparallelopen && taskID ==-1) {
            parallel_access_id = H5Pcreate (H5P_FILE_ACCESS);
            if (parallel_access_id < 0) io_error("Parallel access creation failed");
            herr_t ret = H5Pset_fapl_mpio(parallel_access_id, comm, info);
            if (ret < 0) io_error("Parallel access failed");
            // create the file collectively
            file_id = H5Fcreate(filename.c_str(), flag, H5P_DEFAULT, parallel_access_id);
            if (file_id < 0) io_error(string("Failed to create output file: ")+filename);
            ret = H5Pclose(parallel_access_id);
            if (ret < 0) io_error("Parallel release failed");
            parallel_access_id = -1;
        }
        else {
            if (taskID <0 || taskID > NProcsWrite) io_error(string("MPI Task ID asked to create file out of range. Task ID is ")+to_string(taskID));
            if (ThisWriteTask == taskID) {
                file_id = H5Fcreate(filename.c_str(), flag, H5P_DEFAULT, H5P_DEFAULT);
                if (file_id < 0) io_error(string("Failed to create output file: ")+filename);
                parallel_access_id = -1;
            }
            else {
                parallel_access_id = -2;
            }
            MPI_Barrier(comm);
        }
#else
        file_id = H5Fcreate(filename.c_str(), flag, H5P_DEFAULT, H5P_DEFAULT);
        if(file_id < 0)io_error(string("Failed to create output file: ")+filename);
#endif

    }

    void append(std::string filename, hid_t flag = H5F_ACC_RDWR,
        int taskID = -1, bool iparallelopen = true)
    {
        if(file_id >= 0)io_error("Attempted to open and append to file when already open!");
#ifdef USEPARALLELHDF
        MPI_Comm comm = mpi_comm_write;
        MPI_Info info = MPI_INFO_NULL;
        if (iparallelopen && taskID ==-1) {
            parallel_access_id = H5Pcreate (H5P_FILE_ACCESS);
            if (parallel_access_id < 0) io_error("Parallel access creation failed");
            herr_t ret = H5Pset_fapl_mpio(parallel_access_id, comm, info);
            if (ret < 0) io_error("Parallel access failed");
            // create the file collectively
            file_id = H5Fopen(filename.c_str(), flag, parallel_access_id);
            if (file_id < 0) io_error(string("Failed to create output file: ")+filename);
            ret = H5Pclose(parallel_access_id);
            if (ret < 0) io_error("Parallel release failed");
            parallel_access_id = -1;
        }
        else {
            if (taskID <0 || taskID > NProcsWrite) io_error(string("MPI Task ID asked to create file out of range. Task ID is ")+to_string(taskID));
            if (ThisWriteTask == taskID) {
                file_id = H5Fopen(filename.c_str(),flag, H5P_DEFAULT);
                if (file_id < 0) io_error(string("Failed to create output file: ")+filename);
                parallel_access_id = -1;
            }
            else {
                parallel_access_id = -2;
            }
            MPI_Barrier(comm);
        }
#else
        file_id = H5Fopen(filename.c_str(), flag, H5P_DEFAULT);
        if (file_id < 0) io_error(string("Failed to create output file: ")+filename);
#endif
    }

    // Close the file
    void close()
    {
#ifdef USEPARALLELHDF
        if(file_id < 0 && parallel_access_id == -1) io_error("Attempted to close file which is not open!");
        if (parallel_access_id == -1) H5Fclose(file_id);
#else
        if(file_id < 0) io_error("Attempted to close file which is not open!");
        H5Fclose(file_id);
#endif
        file_id = -1;
#ifdef USEPARALLELHDF
        parallel_access_id = -1;
#endif
    }

    hid_t create_group(string groupname) {
        hid_t group_id = H5Gcreate(file_id, groupname.c_str(),
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        return group_id;
    }
    herr_t close_group(hid_t gid) {
        herr_t status = H5Gclose(gid);
        return status;
    }

  	// Destructor closes the file if it's open
    ~H5OutputFile()
    {
        if(file_id >= 0) close();
    }

    /// Write a new 1D dataset. Data type of the new dataset is taken to be the type of
    /// the input data if not explicitly specified with the filetype_id parameter.
    template <typename T> void write_dataset(std::string name, hsize_t len, T *data,
       hid_t memtype_id = -1, hid_t filetype_id=-1, bool flag_parallel = true, bool flag_hyperslab = true, bool flag_collective = true)
    {
        int rank = 1;
      	hsize_t dims[1] = {len};
        if (memtype_id == -1) memtype_id = hdf5_type(T{});
      	write_dataset_nd(name, rank, dims, data, memtype_id, filetype_id, flag_parallel, flag_hyperslab, flag_collective);
    }
    void write_dataset(string name, hsize_t len, string data, bool flag_parallel = true, bool flag_collective = true)
    {
#ifdef USEPARALLELHDF
        MPI_Comm comm = mpi_comm_write;
        MPI_Info info = MPI_INFO_NULL;
#endif
        int rank = 1;
      	hsize_t dims[1] = {len};

        hid_t memtype_id, filetype_id, dspace_id, dset_id, xfer_plist;
        herr_t status, ret;
        memtype_id = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (memtype_id, data.size());
        filetype_id = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (filetype_id, data.size());

        // Create the dataspace
        dspace_id = H5Screate_simple(rank, dims, NULL);

        // Create the dataset
        dset_id = H5Dcreate(file_id, name.c_str(), filetype_id, dspace_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#ifdef USEPARALLELHDF
        if (flag_parallel) {
            // set up the collective transfer properties list
            xfer_plist = H5Pcreate(H5P_DATASET_XFER);
            if (xfer_plist < 0) io_error(string("Failed to set up parallel transfer: ")+name);
            if (flag_collective) ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
            else ret = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
            if (ret < 0) io_error(string("Failed to set up parallel transfer: ")+name);
            // the result of above should be that all processors write to the same
            // point of the hdf file.
        }
#endif
        // Write the data
        if(H5Dwrite(dset_id, memtype_id, dspace_id, H5S_ALL, H5P_DEFAULT, data.c_str()) < 0)
        io_error(string("Failed to write dataset: ")+name);

        // Clean up (note that dtype_id is NOT a new object so don't need to close it)
#ifdef USEPARALLELHDF
        if (flag_parallel) H5Pclose(xfer_plist);
#endif
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
    }
    void write_dataset(string name, hsize_t len, void *data,
       hid_t memtype_id=-1, hid_t filetype_id=-1, bool flag_parallel = true, bool flag_first_dim_parallel = true, bool flag_hyperslab = true, bool flag_collective = true)
    {
        int rank = 1;
      	hsize_t dims[1] = {len};
        if (memtype_id == -1) {
            throw std::runtime_error("Write data set called with void pointer but no type info passed.");
        }
      	write_dataset_nd(name, rank, dims, data, memtype_id, filetype_id, flag_parallel, flag_first_dim_parallel, flag_hyperslab, flag_collective);
    }


    /// Write a multidimensional dataset. Data type of the new dataset is taken to be the type of
    /// the input data if not explicitly specified with the filetype_id parameter.
    template <typename T> void write_dataset_nd(std::string name, int rank, hsize_t *dims, T *data,
        hid_t memtype_id = -1, hid_t filetype_id = -1,
        bool flag_parallel = true, bool flag_first_dim_parallel = true,
        bool flag_hyperslab = true, bool flag_collective = true)
    {
#ifdef USEPARALLELHDF
        MPI_Comm comm = mpi_comm_write;
        MPI_Info info = MPI_INFO_NULL;
#endif
        hid_t dspace_id, dset_id, prop_id, memspace_id, ret;
        vector<hsize_t> chunks(rank);

        // Get HDF5 data type of the array in memory
        if (memtype_id == -1) memtype_id = hdf5_type(T{});

        // Determine type of the dataset to create
        if(filetype_id < 0) filetype_id = memtype_id;

#ifdef USEPARALLELHDF
        vector<unsigned long long> mpi_hdf_dims(rank*NProcsWrite), mpi_hdf_dims_tot(rank), dims_single(rank), dims_offset(rank);
        if (flag_parallel) {
            //if parallel hdf5 get the full extent of the data
            //this bit of code communicating information can probably be done elsewhere
            //minimize number of mpi communications
            for (auto i=0;i<rank;i++) dims_single[i]=dims[i];
            MPI_Allgather(dims_single.data(), rank, MPI_UNSIGNED_LONG_LONG, mpi_hdf_dims.data(), rank, MPI_UNSIGNED_LONG_LONG, comm);
            MPI_Allreduce(dims_single.data(), mpi_hdf_dims_tot.data(), rank, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
            for (auto i=0;i<rank;i++) {
                dims_offset[i] = 0;
                if (flag_first_dim_parallel && i > 0) continue;
                for (auto j=1;j<=ThisWriteTask;j++) {
                    dims_offset[i] += mpi_hdf_dims[i*NProcs+j-1];
                }
            }
            if (flag_first_dim_parallel && rank > 1) {
                for (auto i=1; i<rank;i++) mpi_hdf_dims_tot[i] = dims[i];
            }
        }
#endif

        // Determine if going to compress data in chunks
        // Only chunk non-zero size datasets
        int nonzero_size = 1;
        for(int i=0; i<rank; i++)
        {
#ifdef USEPARALLELHDF
            if (flag_parallel) {
                if(mpi_hdf_dims_tot[i]==0) nonzero_size = 0;
            }
            else {
                if(dims[i]==0) nonzero_size = 0;
            }
#else
            if(dims[i]==0) nonzero_size = 0;
#endif
        }
        // Only chunk datasets where we would have >1 chunk
        int large_dataset = 0;
        for(int i=0; i<rank; i++)
        {
#ifdef USEPARALLELHDF
            if (flag_parallel) {
                if(mpi_hdf_dims_tot[i] > HDFOUTPUTCHUNKSIZE) large_dataset = 1;
            }
            else {
                if(dims[i] > HDFOUTPUTCHUNKSIZE) large_dataset = 1;
            }
#else
            if(dims[i] > HDFOUTPUTCHUNKSIZE) large_dataset = 1;
#endif
        }
        if(nonzero_size && large_dataset)
        {
#ifdef USEPARALLELHDF
            if (flag_parallel) {
                for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, mpi_hdf_dims_tot[i]);
            }
            else {
                for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, dims[i]);
            }
#else
            for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, dims[i]);
#endif
        }

        // Create the dataspace
#ifdef USEPARALLELHDF
        if (flag_parallel) {
            //then all threads create the same simple data space
            //so the meta information is the same
            if (flag_hyperslab) {
                //allocate the space spanning the file
                dspace_id = H5Screate_simple(rank, mpi_hdf_dims_tot.data(), NULL);
                //allocate the memory space
                //allocate the memory space
                memspace_id = H5Screate_simple(rank, dims, NULL);
            }
            else {
                dspace_id = H5Screate_simple(rank, dims, NULL);
                memspace_id = dspace_id;
            }
        }
        else {
            dspace_id = H5Screate_simple(rank, dims, NULL);
            memspace_id = dspace_id;
        }
#else
        dspace_id = H5Screate_simple(rank, dims, NULL);
        memspace_id = dspace_id;
#endif

        // Dataset creation properties
        prop_id = H5P_DEFAULT;
#ifdef USEHDFCOMPRESSOIN
        // this defines compression
        if(nonzero_size && large_dataset)
        {
            prop_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_layout(prop_id, H5D_CHUNKED);
            H5Pset_chunk(prop_id, rank, chunks.data());
            H5Pset_deflate(prop_id, HDFDEFLATE);
        }
#endif
        // Create the dataset
        dset_id = H5Dcreate(file_id, name.c_str(), filetype_id, dspace_id,
            H5P_DEFAULT, prop_id, H5P_DEFAULT);
        if(dset_id < 0) io_error(string("Failed to create dataset: ")+name);
        H5Pclose(prop_id);

        prop_id = H5P_DEFAULT;
#ifdef USEPARALLELHDF
        if (flag_parallel) {
            // set up the collective transfer properties list
            prop_id = H5Pcreate(H5P_DATASET_XFER);
            //if all tasks are participating in the writes
            if (flag_collective) ret = H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
            else ret = H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_INDEPENDENT);
            if (flag_hyperslab) {
                H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, dims_offset.data(), NULL, dims, NULL);
                if (dims[0] == 0) {
                    H5Sselect_none(dspace_id);
                    H5Sselect_none(memspace_id);
                }

            }
            if (mpi_hdf_dims_tot[0] > 0) {
                // Write the data
                ret = H5Dwrite(dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
                if (ret < 0) io_error(string("Failed to write dataset: ")+name);
            }
        }
        else if (dims[0] > 0)
        {
            // Write the data
            ret = H5Dwrite(dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
            if (ret < 0) io_error(string("Failed to write dataset: ")+name);
        }

#else
        // Write the data
        if (dims[0] > 0) {
            ret = H5Dwrite(dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
            if (ret < 0) io_error(string("Failed to write dataset: ")+name);
        }
#endif

        // Clean up (note that dtype_id is NOT a new object so don't need to close it)
        H5Pclose(prop_id);
#ifdef USEPARALLELHDF
        if (flag_hyperslab && flag_parallel) H5Sclose(memspace_id);
#endif
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
    }
    void write_dataset_nd(std::string name, int rank, hsize_t *dims, void *data,
        hid_t memtype_id = -1, hid_t filetype_id=-1,
        bool flag_parallel = true, bool flag_first_dim_parallel = true,
        bool flag_hyperslab = true, bool flag_collective = true)
    {
#ifdef USEPARALLELHDF
        MPI_Comm comm = mpi_comm_write;
        MPI_Info info = MPI_INFO_NULL;
#endif
        hid_t dspace_id, dset_id, prop_id, memspace_id, ret;
        vector<hsize_t> chunks(rank);
        // Get HDF5 data type of the array in memory
        if (memtype_id == -1) {
            throw std::runtime_error("Write data set called with void pointer but no type info passed.");
        }
        // Determine type of the dataset to create
        if(filetype_id < 0) filetype_id = memtype_id;

#ifdef USEPARALLELHDF
        vector<unsigned long long> mpi_hdf_dims(rank*NProcsWrite), mpi_hdf_dims_tot(rank), dims_single(rank), dims_offset(rank);
        //if parallel hdf5 get the full extent of the data
        //this bit of code communicating information can probably be done elsewhere
        //minimize number of mpi communications
        if (flag_parallel) {
            //if parallel hdf5 get the full extent of the data
            //this bit of code communicating information can probably be done elsewhere
            //minimize number of mpi communications
            for (auto i=0;i<rank;i++) dims_single[i]=dims[i];
            MPI_Allgather(dims_single.data(), rank, MPI_UNSIGNED_LONG_LONG, mpi_hdf_dims.data(), rank, MPI_UNSIGNED_LONG_LONG, comm);
            MPI_Allreduce(dims_single.data(), mpi_hdf_dims_tot.data(), rank, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
            for (auto i=0;i<rank;i++) {
                dims_offset[i] = 0;
                if (flag_first_dim_parallel && i > 0) continue;
                for (auto j=1;j<=ThisWriteTask;j++) {
                    dims_offset[i] += mpi_hdf_dims[i*NProcs+j-1];
                }
            }
            if (flag_first_dim_parallel && rank > 1) {
                for (auto i=1; i<rank;i++) mpi_hdf_dims_tot[i] = dims[i];
            }
        }
#endif

        // Determine if going to compress data in chunks
        // Only chunk non-zero size datasets
        int nonzero_size = 1;
        for(int i=0; i<rank; i++)
        {
#ifdef USEPARALLELHDF
            if (flag_parallel) {
                if(mpi_hdf_dims_tot[i]==0) nonzero_size = 0;
            }
            else {
                if(dims[i]==0) nonzero_size = 0;
            }
#else
            if(dims[i]==0) nonzero_size = 0;
#endif
        }
        // Only chunk datasets where we would have >1 chunk
        int large_dataset = 0;
        for(int i=0; i<rank; i++)
        {
#ifdef USEPARALLELHDF
            if (flag_parallel) {
                if(mpi_hdf_dims_tot[i] > HDFOUTPUTCHUNKSIZE) large_dataset = 1;
            }
            else {
                if(dims[i] > HDFOUTPUTCHUNKSIZE) large_dataset = 1;
            }
#else
            if(dims[i] > HDFOUTPUTCHUNKSIZE) large_dataset = 1;
#endif
        }
        if(nonzero_size && large_dataset)
        {
#ifdef USEPARALLELHDF
            if (flag_parallel) {
                for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, mpi_hdf_dims_tot[i]);
            }
            else {
                for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, dims[i]);
            }
#else
            for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, dims[i]);
#endif
        }

        // Create the dataspace
#ifdef USEPARALLELHDF
        if (flag_parallel) {
            //then all threads create the same simple data space
            //so the meta information is the same
            if (flag_hyperslab) {
                //allocate the space spanning the file
                dspace_id = H5Screate_simple(rank, mpi_hdf_dims_tot.data(), NULL);
                //allocate the memory space
                memspace_id = H5Screate_simple(rank, dims, NULL);
            }
            else {
                dspace_id = H5Screate_simple(rank, dims, NULL);
                memspace_id = dspace_id;
            }
        }
        else {
            dspace_id = H5Screate_simple(rank, dims, NULL);
            memspace_id = dspace_id;
        }
#else
        dspace_id = H5Screate_simple(rank, dims, NULL);
        memspace_id = dspace_id;
#endif

        // Dataset creation properties
        prop_id = H5P_DEFAULT;
#ifdef USEHDFCOMPRESSOIN
        // this defines compression
        if(nonzero_size && large_dataset)
        {
            prop_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_layout(prop_id, H5D_CHUNKED);
            H5Pset_chunk(prop_id, rank, chunks.data());
            H5Pset_deflate(prop_id, HDFDEFLATE);
        }
#endif

        // Create the dataset
        dset_id = H5Dcreate(file_id, name.c_str(), filetype_id, dspace_id,
            H5P_DEFAULT, prop_id, H5P_DEFAULT);
        if(dset_id < 0) io_error(string("Failed to create dataset: ")+name);
        H5Pclose(prop_id);

        prop_id = H5P_DEFAULT;

#ifdef USEPARALLELHDF
        if (flag_parallel) {
            // set up the collective transfer properties list
            prop_id = H5Pcreate(H5P_DATASET_XFER);
            //if all tasks are participating in the writes
            if (flag_collective) ret = H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_COLLECTIVE);
            else ret = H5Pset_dxpl_mpio(prop_id, H5FD_MPIO_INDEPENDENT);
            if (flag_hyperslab) {
                H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, dims_offset.data(), NULL, dims, NULL);
                if (dims[0] == 0) {
                    H5Sselect_none(dspace_id);
                    H5Sselect_none(memspace_id);
                }

            }
            if (mpi_hdf_dims_tot[0] > 0) {
                // Write the data
                ret = H5Dwrite(dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
                if (ret < 0) io_error(string("Failed to write dataset: ")+name);
            }
        }
        else if (dims[0] > 0)
        {
            // Write the data
            ret = H5Dwrite(dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
            if (ret < 0) io_error(string("Failed to write dataset: ")+name);
        }

#else
        // Write the data
        if (dims[0] > 0) {
            ret = H5Dwrite(dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
            if (ret < 0) io_error(string("Failed to write dataset: ")+name);
        }
#endif

        // Clean up (note that dtype_id is NOT a new object so don't need to close it)
        H5Pclose(prop_id);
#ifdef USEPARALLELHDF
        if (flag_hyperslab && flag_parallel) H5Sclose(memspace_id);
#endif
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
    }

    /// write an attribute
    template <typename T> void write_attribute(std::string parent, std::string name, T data)
    {
        // Get HDF5 data type of the value to write
        hid_t dtype_id = hdf5_type(data);

        // Open the parent object
        hid_t parent_id = H5Oopen(file_id, parent.c_str(), H5P_DEFAULT);
        if(parent_id < 0)io_error(string("Unable to open object to write attribute: ")+name);

        // Create dataspace
        hid_t dspace_id = H5Screate(H5S_SCALAR);

        // Create attribute
        hid_t attr_id = H5Acreate(parent_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
        if(attr_id < 0)io_error(string("Unable to create attribute ")+name+string(" on object ")+parent);

        // Write the attribute
        if(H5Awrite(attr_id, dtype_id, &data) < 0)
        io_error(string("Unable to write attribute ")+name+string(" on object ")+parent);

        // Clean up
        H5Aclose(attr_id);
        H5Sclose(dspace_id);
        H5Oclose(parent_id);
    }

    void write_attribute(string parent, string name, string data)
    {
        // Get HDF5 data type of the value to write
        hid_t dtype_id = H5Tcopy(H5T_C_S1);
        if (data.size() == 0) data=" ";
        H5Tset_size(dtype_id, data.size());
        H5Tset_strpad(dtype_id, H5T_STR_NULLTERM);

        // Open the parent object
        hid_t parent_id = H5Oopen(file_id, parent.c_str(), H5P_DEFAULT);
        if(parent_id < 0)io_error(string("Unable to open object to write attribute: ")+name);

        // Create dataspace
        hid_t dspace_id = H5Screate(H5S_SCALAR);

        // Create attribute
        hid_t attr_id = H5Acreate(parent_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
        if(attr_id < 0)io_error(string("Unable to create attribute ")+name+string(" on object ")+parent);

        // Write the attribute
        if(H5Awrite(attr_id, dtype_id, data.c_str()) < 0)
        io_error(string("Unable to write attribute ")+name+string(" on object ")+parent);

        // Clean up
        H5Aclose(attr_id);
        H5Sclose(dspace_id);
        H5Oclose(parent_id);
    }
};


///This structures stores the strings defining the groups of data in the hdf input. NOTE: HERE I show the strings for Illustris format
struct HDF_Group_Names {
    //define the strings associated with the types of structures contained in the hdf file.
    string Header_name;
    string GASpart_name;
    string DMpart_name;
    string EXTRADMpart_name;
    string EXTRApart_name;
    string TRACERpart_name;
    string STARpart_name;
    string BHpart_name;
    string part_names[NHDFTYPE];
    string names[NHDFTYPE+1];

    ///constructor
    HDF_Group_Names(int hdfnametype=HDFEAGLENAMES){
        switch (hdfnametype) {
          case HDFSWIFTEAGLENAMES:
            Header_name=string("Header");
            GASpart_name=string("PartType0");
            DMpart_name=string("PartType1");
            EXTRADMpart_name=string("PartType2");
            EXTRApart_name=string("PartType2");
            TRACERpart_name=string("PartType3");
            STARpart_name=string("PartType4");
            BHpart_name=string("PartType5");
          break;

          default:
            Header_name=string("Header");
            GASpart_name=string("PartType0");
            DMpart_name=string("PartType1");
            EXTRADMpart_name=string("PartType2");
            EXTRApart_name=string("PartType2");
            TRACERpart_name=string("PartType3");
            STARpart_name=string("PartType4");
            BHpart_name=string("PartType5");
          break;
        }

        part_names[0]=GASpart_name;
        part_names[1]=DMpart_name;
#ifdef HIGHRES
        part_names[2]=EXTRADMpart_name;
#else
        part_names[2]=EXTRApart_name;
#endif
        part_names[3]=TRACERpart_name;
        part_names[4]=STARpart_name;
        part_names[5]=BHpart_name;

        names[0]=Header_name;
        names[1]=GASpart_name;
        names[2]=DMpart_name;
#ifdef HIGHRES
        names[3]=EXTRADMpart_name;
#else
        names[3]=EXTRApart_name;
#endif
        names[4]=TRACERpart_name;
        names[5]=STARpart_name;
        names[6]=BHpart_name;
    }
};

///data stored in the header group structure in the HDF format
struct HDF_Header {

    double      BoxSize;
    unsigned long long npart[NHDFTYPE];
    unsigned int npartTotal[NHDFTYPE];
    unsigned int npartTotalHW[NHDFTYPE];
    double      mass[NHDFTYPE];
    double      Omega0, OmegaLambda, HubbleParam;
    double      redshift, time;
    int         iscosmological;
    int         num_files;

    string names[HDFHEADNINFO];
    const static int IBoxSize  =0;
    const static int IMass     =1;
    const static int INuminFile=2;
    const static int INumTot   =3;
    const static int INumTotHW =4;
    const static int IOmega0   =5;
    const static int IOmegaL   =6;
    const static int IRedshift =7;
    const static int ITime     =8;
    const static int INumFiles =9;
    const static int IHubbleParam =10;
    const static int IIsCosmological =11;

    ///constructor
    HDF_Header(int hdfnametype=HDFEAGLENAMES) {
        int itemp=0;
        switch (hdfnametype) {
          case HDFSWIFTEAGLENAMES:
            names[itemp++]=string("Header/BoxSize");
            names[itemp++]=string("Header/MassTable");
            names[itemp++]=string("Header/NumPart_ThisFile");
            names[itemp++]=string("Header/NumPart_Total");
            names[itemp++]=string("Header/NumPart_Total_HighWord");
            names[itemp++]=string("Cosmology/Omega_m");
            names[itemp++]=string("Cosmology/Omega_lambda");
            names[itemp++]=string("Header/Redshift");
            names[itemp++]=string("Header/Time");
            names[itemp++]=string("Header/NumFilesPerSnapshot");
            names[itemp++]=string("Cosmology/h");
            names[itemp++]=string("Cosmology/Cosmological run");
            break;
          case HDFOLDSWIFTEAGLENAMES:
            names[itemp++]=string("Header/BoxSize");
            names[itemp++]=string("Header/MassTable");
            names[itemp++]=string("Header/NumPart_ThisFile");
            names[itemp++]=string("Header/NumPart_Total");
            names[itemp++]=string("Header/NumPart_Total_HighWord");
            names[itemp++]=string("Cosmology/Omega_m");
            names[itemp++]=string("Cosmology/Omega_lambda");
            names[itemp++]=string("Header/Redshift");
            names[itemp++]=string("Header/Time");
            names[itemp++]=string("Header/NumFilesPerSnapshot");
            names[itemp++]=string("Cosmology/h");
            names[itemp++]=string("Cosmology/Cosmological run");
            break;

          default:
            names[itemp++]=string("Header/BoxSize");
            names[itemp++]=string("Header/MassTable");
            names[itemp++]=string("Header/NumPart_ThisFile");
            names[itemp++]=string("Header/NumPart_Total");
            names[itemp++]=string("Header/NumPart_Total_HighWord");
            names[itemp++]=string("Header/Omega0");
            names[itemp++]=string("Header/OmegaLambda");
            names[itemp++]=string("Header/Redshift");
            names[itemp++]=string("Header/Time");
            names[itemp++]=string("Header/NumFilesPerSnapshot");
            names[itemp++]=string("Header/HubbleParam");
            break;
        }
    }
};

struct HDF_Part_Info {
    string names[HDFMAXPINFO];
    int ptype;
    int nentries;
    //store where properties are located
    int propindex[100];

    //the HDF naming convenction for the data blocks. By default assumes ILLUSTRIS nameing convention
    //for simplicity, all particles have basic properties listed first, x,v,ids,mass in this order
    HDF_Part_Info(int PTYPE, int hdfnametype=HDFEAGLENAMES) {
        ptype=PTYPE;
        int itemp=0;
        //gas
        if (ptype==HDFGASTYPE) {

            // Positions
            names[itemp++]=string("Coordinates");

            // Velocities
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");

            // IDs
            names[itemp++]=string("ParticleIDs");

            // Masses
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Mass");
            else names[itemp++]=string("Masses");

            // Density
            if(hdfnametype==HDFSWIFTEAGLENAMES) names[itemp++]=string("Densities");
            else names[itemp++]=string("Density");

            // Internal energies
            if(hdfnametype==HDFSWIFTEAGLENAMES) names[itemp++]=string("InternalEnergies");
            else names[itemp++]=string("InternalEnergy");

            // SFR
            if(hdfnametype==HDFSWIFTEAGLENAMES) names[itemp++]=string("StarFormationRates");
            else if(hdfnametype==HDFOLDSWIFTEAGLENAMES) names[itemp++]=string("SFR");
            else names[itemp++]=string("StarFormationRate");

            //Metallicity. Note always place at position 7 in naming array
            if (hdfnametype==HDFILLUSTISNAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=string("GFM_Metallicity");
                names[itemp++]=string("ElectronAbundance");
                names[itemp++]=string("NeutralHydrogenAbundance");
                names[itemp++]=string("Volume");
                names[itemp++]=string("SmoothingLength");
                names[itemp++]=string("Potential");
                names[itemp++]=string("SubfindDensity");
                names[itemp++]=string("SubfindHsml");
                names[itemp++]=string("SubfindVelDisp");
                names[itemp++]=string("GFM_AGNRadiation");
                names[itemp++]=string("GFM_CoolingRate");
                names[itemp++]=string("GFM_WindDMVelDisp");
                names[itemp++]=string("NumTracers");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                propindex[HDFGASIMETAL]=itemp;
                //names[itemp++]=string("Metallicity");//11 metals stored in this data set
                names[itemp++]=string("Metallicity_00");//only grab the first of the 11, which is total
                names[itemp++]=string("ElectronAbundance");
                names[itemp++]=string("FractionH2");
                names[itemp++]=string("GrackleHI");
                names[itemp++]=string("GrackleHII");
                names[itemp++]=string("GrackleHM");
                names[itemp++]=string("GrackleHeI");
                names[itemp++]=string("GrackleHeII");
                names[itemp++]=string("GrackleHeIII");
                names[itemp++]=string("NWindLaunches");
                names[itemp++]=string("NeutralHydrogenAbundance");
                names[itemp++]=string("ParticleChildIDsNumber");
                names[itemp++]=string("ParticleIDGenerationNumber");
                names[itemp++]=string("Potential");
                names[itemp++]=string("Sigma");
                names[itemp++]=string("SmoothingLength");
            }
            else if (hdfnametype==HDFSIMBANAMES || hdfnametype==HDFMUFASANAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=string("Metallicity");//11 metals stored in this data set
                names[itemp++]=string("ElectronAbundance");
                names[itemp++]=string("FractionH2");
                names[itemp++]=string("GrackleHI");
                names[itemp++]=string("GrackleHII");
                names[itemp++]=string("GrackleHM");
                names[itemp++]=string("GrackleHeI");
                names[itemp++]=string("GrackleHeII");
                names[itemp++]=string("GrackleHeIII");
                names[itemp++]=string("NWindLaunches");
                names[itemp++]=string("NeutralHydrogenAbundance");
                names[itemp++]=string("Potential");
                names[itemp++]=string("Sigma");
                names[itemp++]=string("SmoothingLength");
                names[itemp++]=string("DelayTime");
                names[itemp++]=string("Dust_Masses");
                names[itemp++]=string("Dust_Metallicity");//11 metals stored in this data set
            }
            else if (hdfnametype==HDFEAGLENAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=string("Metallicity");
            }
            else if(hdfnametype==HDFSWIFTEAGLENAMES) {
              propindex[HDFGASIMETAL]=itemp;
              names[itemp++]=string("MetalMassFractions");
            }
            else if(hdfnametype==HDFOLDSWIFTEAGLENAMES) {
              propindex[HDFGASIMETAL]=itemp;
              names[itemp++]=string("Metallicity");
            }

        }
        //dark matter
        if (ptype==HDFDMTYPE) {

            // Positions
            names[itemp++]=string("Coordinates");

            // Velocities
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");

            // IDs
            names[itemp++]=string("ParticleIDs");

            // Masses
            if (hdfnametype==HDFSWIFTEAGLENAMES || hdfnametype==HDFSIMBANAMES ||
                hdfnametype==HDFMUFASANAMES || hdfnametype==HDFOLDSWIFTEAGLENAMES) {
                names[itemp++]=string("Masses");
            }

            // Potential
            if (hdfnametype==HDFSIMBANAMES||hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=string("Potential");
            }
            else if (hdfnametype==HDFILLUSTISNAMES) {
                names[itemp++]=string("Potential");
            }

            // Subfind properties
            if (hdfnametype==HDFILLUSTISNAMES) {
                names[itemp++]=string("SubfindDensity");
                names[itemp++]=string("SubfindHsml");
                names[itemp++]=string("SubfindVelDisp");
            }
        }

        //also dark matter particles
        if (ptype==HDFDM1TYPE ||ptype==HDFDM2TYPE) {

            // Positions
            names[itemp++]=string("Coordinates");

            // Velocities
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");

            // IDs
            names[itemp++]=string("ParticleIDs");

            // Masses
            names[itemp++]=string("Masses");

            // Potential
            if (hdfnametype==HDFSIMBANAMES||hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=string("Potential");
            }
        }
        if (ptype==HDFTRACERTYPE) {
            names[itemp++]=string("FluidQuantities");
            names[itemp++]=string("ParentID");
            names[itemp++]=string("TracerID");
        }
        if (ptype==HDFSTARTYPE) {

            // Positions
            names[itemp++]=string("Coordinates");

            // Velocities
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");

            // IDs
            names[itemp++]=string("ParticleIDs");

            // Masses
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Mass");
            else names[itemp++]=string("Masses");

            //for stars assume star formation and metallicy are position 4, 5 in name array
            if (hdfnametype==HDFILLUSTISNAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=string("GFM_StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=string("GFM_Metallicity");
                names[itemp++]=string("Potential");
                names[itemp++]=string("SubfindDensity");
                names[itemp++]=string("SubfindHsml");
                names[itemp++]=string("SubfindVelDisp");
                names[itemp++]=string("GFM_InitialMass");
                names[itemp++]=string("GFM_StellarPhotometrics");
                names[itemp++]=string("NumTracers");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=string("StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                //names[itemp++]=string("Metallicity");//11 metals stored in this data set
                names[itemp++]=string("Metallicity_00");//only grab the first of the 11, which is total
                names[itemp++]=string("AGS-Softening Dataset");
                names[itemp++]=string("ParticleChildIDsNumber");
                names[itemp++]=string("ParticleIDGenerationNumber");
                names[itemp++]=string("Potential");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=string("StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=string("Metallicity");
                names[itemp++]=string("Potential");
                names[itemp++]=string("Dust_Masses");
                names[itemp++]=string("Dust_Metallicity");//11 metals stored in this data set
            }
            else if (hdfnametype==HDFEAGLENAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=string("StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=string("Metallicity");
            }
            else if (hdfnametype==HDFSWIFTEAGLENAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=string("BirthScaleFactors");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=string("MetalMassFractions");
            }
        }
        if (ptype==HDFBHTYPE) {

            // Positions
            names[itemp++]=string("Coordinates");

            // Velocities
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");

            // IDs
            names[itemp++]=string("ParticleIDs");

            // Masses
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Mass");
            if(hdfnametype==HDFSWIFTEAGLENAMES) names[itemp++]=string("DynamicalMasses");
            else names[itemp++]=string("Masses");

            if (hdfnametype==HDFILLUSTISNAMES) {
                names[itemp++]=string("HostHaloMass");
                names[itemp++]=string("Potential");
                names[itemp++]=string("SubfindDensity");
                names[itemp++]=string("SubfindHsml");
                names[itemp++]=string("SubfindVelDisp");
                names[itemp++]=string("BH_CumEgyInjection_QM");
                names[itemp++]=string("BH_CumMassGrowth_QM");
                names[itemp++]=string("BH_Density");
                names[itemp++]=string("BH_Hsml");
                names[itemp++]=string("BH_Mass");
                names[itemp++]=string("BH_Mass_bubbles");
                names[itemp++]=string("BH_Mass_ini");
                names[itemp++]=string("BH_Mdot");
                names[itemp++]=string("BH_Pressure");
                names[itemp++]=string("BH_Progs");
                names[itemp++]=string("BH_U");
                names[itemp++]=string("NumTracers");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                names[itemp++]=string("AGS-Softening");
                names[itemp++]=string("Potential");
            }
            else if (hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=string("StellarFormationTime");
                names[itemp++]=string("BH_AccretionLength");
                names[itemp++]=string("BH_Mass");
                names[itemp++]=string("BH_Mass_AlphaDisk");
                names[itemp++]=string("BH_Mass_Mdot");
                names[itemp++]=string("BH_NProgs");
                names[itemp++]=string("Potential");
            }
            else if (hdfnametype==HDFEAGLENAMES) {
                //names[itemp++]=string("StellarFormationTime");
                //names[itemp++]=string("Metallicity");
            }
            else if (hdfnametype==HDFSWIFTEAGLENAMES) {
                propindex[HDFBHIAGE]=itemp;
                names[itemp++]=string("FormationScaleFactors");
                propindex[HDFBHIMETAL]=itemp;
                names[itemp++]=string("MetalMasses");
                propindex[HDFBHIMDOT]=itemp;
                names[itemp++]=string("AccretionRates");

                names[itemp++]=string("SubgridMasses");
                names[itemp++]=string("ElementMasses");
                names[itemp++]=string("MetalMassFromSNIa");
                names[itemp++]=string("MetalMassFromSNII");
                names[itemp++]=string("MetalMassFromAGB");
                names[itemp++]=string("MassesFromSNIa");
                names[itemp++]=string("MassesFromSNII");
                names[itemp++]=string("MassesFromAGB");
                names[itemp++]=string("IronMassFromSNIa");
                names[itemp++]=string("GasDensities");
                names[itemp++]=string("GasSoundSpeeds");
                names[itemp++]=string("EnergyReservoirs");
                names[itemp++]=string("TotalAccretedMasses");
            }
        }
        nentries=itemp;
    }
};
//@}

/// \name Set the particle types to be load in from the HDF file
//@{
inline void HDFSetUsedParticleTypes(Options &opt, int &nusetypes, int &nbusetypes, int usetypes[])
{
    nusetypes=0;
    if (opt.partsearchtype==PSTALL) {
        nusetypes=0;
        if (opt.iusegasparticles) usetypes[nusetypes++]=HDFGASTYPE;
        if (opt.iusedmparticles) usetypes[nusetypes++]=HDFDMTYPE;
        if (opt.iuseextradarkparticles) {
            usetypes[nusetypes++]=HDFDM1TYPE;
            if (opt.ihdfnameconvention!=HDFSWIFTEAGLENAMES)
            {
                usetypes[nusetypes++]=HDFDM2TYPE;
            }
        }
        if (opt.iusestarparticles) usetypes[nusetypes++]=HDFSTARTYPE;
        if (opt.iusesinkparticles) usetypes[nusetypes++]=HDFBHTYPE;
        if (opt.iusewindparticles) usetypes[nusetypes++]=HDFWINDTYPE;
        if (opt.iusetracerparticles) usetypes[nusetypes++]=HDFTRACERTYPE;
    }
    else if (opt.partsearchtype==PSTDARK) {
        nusetypes=1;usetypes[0]=HDFDMTYPE;
        if (opt.iuseextradarkparticles) {
            usetypes[nusetypes++]=HDFDM1TYPE;
            if (opt.ihdfnameconvention!=HDFSWIFTEAGLENAMES)
            {
                usetypes[nusetypes++]=HDFDM2TYPE;
            }
        }
        if (opt.iBaryonSearch) {
            nbusetypes=1;usetypes[nusetypes+nbusetypes++]=HDFGASTYPE;
            if (opt.iusestarparticles) usetypes[nusetypes+nbusetypes++]=HDFSTARTYPE;
            if (opt.iusesinkparticles) usetypes[nusetypes+nbusetypes++]=HDFBHTYPE;
        }
    }
    else if (opt.partsearchtype==PSTGAS) {nusetypes=1;usetypes[0]=HDFGASTYPE;}
    else if (opt.partsearchtype==PSTSTAR) {nusetypes=1;usetypes[0]=HDFSTARTYPE;}
    else if (opt.partsearchtype==PSTBH) {nusetypes=1;usetypes[0]=HDFBHTYPE;}
}
//@}

/// \name Get the number of particles in the hdf files
//@{
inline Int_t HDF_get_nbodies(char *fname, int ptype, Options &opt)
{
    char buf[2000],buf1[2000],buf2[2000];
    sprintf(buf1,"%s.0.hdf5",fname);
    sprintf(buf2,"%s.hdf5",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    //H5File Fhdf;
    hid_t Fhdf;
    HDF_Group_Names hdf_gnames;
    //to store the groups, data sets and their associated data spaces
    //Attribute headerattribs;
    hid_t headerattribs;
    HDF_Header hdf_header_info = HDF_Header(opt.ihdfnameconvention);
    //buffers to load data
    string stringbuff, dataname;
    string swift_str = "SWIFT";
    int intbuff[NHDFTYPE];
    long long longbuff[NHDFTYPE];
    unsigned int uintbuff[NHDFTYPE];
    vector<unsigned int> vuintbuff;
    int j,k,ireaderror=0;
    Int_t nbodies=0;
    //DataSpace headerdataspace;
    hid_t headerdataspace;

    //to determine types
    //IntType inttype;
    //StrType stringtype;
    int nusetypes,usetypes[NHDFTYPE],nbusetypes;
    HDFSetUsedParticleTypes(opt,nusetypes,nbusetypes,usetypes);

    //Try block to detect exceptions raised by any of the calls inside it
    //try
    {
        //turn off the auto-printing when failure occurs so that we can
        //handle the errors appropriately
        //Exception::dontPrint();

        //Open the specified file and the specified dataset in the file.
        //Fhdf.openFile(buf, H5F_ACC_RDONLY);
        Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
        cout<<"Loading HDF header info in header group: "<<hdf_gnames.Header_name<<endl;

        if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES) {

            // Check if it is a SWIFT snapshot.
            //headerattribs=get_attribute(Fhdf, "Header/Code");
            //stringtype = headerattribs.getStrType();
            //headerattribs.read(stringtype, stringbuff);
            dataname = string("Header/Code");
            stringbuff = read_attribute<string>(Fhdf, dataname);

            // Read SWIFT parameters
            if(!swift_str.compare(stringbuff)) {
                hdf_header_info.iscosmological=read_attribute<int>(Fhdf, hdf_header_info.names[hdf_header_info.IIsCosmological]);
                if (!hdf_header_info.iscosmological && opt.icosmologicalin) {
                    cout<<"Error: cosmology is turned on in the config file but the snaphot provided is a non-cosmological run."<<endl;
#ifdef USEMPI
                    MPI_Abort(MPI_COMM_WORLD, 8);
#else
                    exit(0);
#endif
                }
                else if (hdf_header_info.iscosmological && !opt.icosmologicalin) {
                    cout<<"Error: cosmology is turned off in the config file but the snaphot provided is a cosmological run."<<endl;
#ifdef USEMPI
                    MPI_Abort(MPI_COMM_WORLD, 8);
#else
                    exit(0);
#endif
                }
            }
            // If the code is not SWIFT
            else {
                cout<<"SWIFT EAGLE HDF5 naming convention chosen in config file but the snapshot was not produced by SWIFT. The string read was: "<<stringbuff<<endl;
#ifdef USEMPI
                MPI_Abort(MPI_COMM_WORLD, 8);
#else
                exit(0);
#endif
            }
        }

        vuintbuff=read_attribute_v<unsigned int>(Fhdf, hdf_header_info.names[hdf_header_info.INumTot]);
        for (j=0;j<NHDFTYPE;j++) hdf_header_info.npartTotal[j]=vuintbuff[j];
        vuintbuff=read_attribute_v<unsigned int>(Fhdf, hdf_header_info.names[hdf_header_info.INumTotHW]);
        for (j=0;j<NHDFTYPE;j++) hdf_header_info.npartTotalHW[j]=vuintbuff[j];

        if (opt.partsearchtype==PSTALL) {
            if (opt.iusestarparticles && hdf_header_info.npartTotalHW[HDFSTARTYPE] == 0 && hdf_header_info.npartTotal[HDFSTARTYPE] == 0)
            {
                cerr<<"Warning: Configured to load star particles but none present"<<endl;
                opt.iusestarparticles=0;
            }
            if (opt.iusesinkparticles && hdf_header_info.npartTotalHW[HDFBHTYPE] == 0 && hdf_header_info.npartTotal[HDFBHTYPE] == 0)
            {
                cerr<<"Warning: Configured to load black hole particles but none present"<<endl;
                opt.iusesinkparticles=0;
            }
            if (opt.iusewindparticles && hdf_header_info.npartTotalHW[HDFWINDTYPE] == 0 && hdf_header_info.npartTotal[HDFWINDTYPE] == 0)
            {
                cerr<<"Warning: Configured to load wind particles but none present"<<endl;
                opt.iusewindparticles=0;
            }
            if (opt.iusetracerparticles && hdf_header_info.npartTotalHW[HDFTRACERTYPE] == 0 && hdf_header_info.npartTotal[HDFTRACERTYPE] == 0)
            {
                cerr<<"Warning: Configured to load tracer particles but none present"<<endl;
                opt.iusetracerparticles=0;
            }
#ifdef HIGHRES
            if (opt.iuseextradarkparticles && hdf_header_info.npartTotalHW[HDFDM1TYPE] == 0 && hdf_header_info.npartTotal[HDFDM1TYPE] == 0
                && hdf_header_info.npartTotalHW[HDFDM2TYPE] == 0 && hdf_header_info.npartTotal[HDFDM2TYPE] == 0)
            {
                cerr<<"Warning: Configured to load extra low res dark matter particles but none present"<<endl;
                opt.iuseextradarkparticles=0;
            }
#endif
        }
    }
    /*
    catch(GroupIException &error)
    {
        HDF5PrintError(error);
        cerr<<"Error in group might suggest config file has the incorrect HDF naming convention. ";
        cerr<<"Check HDF_name_convetion or add new naming convention updating hdfitems.h in the source code. "<<endl;
        Fhdf.close();
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
      Fhdf.close();
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
      Fhdf.close();
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
      Fhdf.close();
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
      Fhdf.close();
#ifdef USEMPI
      MPI_Abort(MPI_COMM_WORLD,8);
#else
      exit(8);
#endif
    }
    // catch failure caused by missing attribute
    catch( invalid_argument error )
    {
      if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) {
        cerr<<"Reading SWIFT EAGLE HDF5 file: "<<error.what()<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD, 8);
#else
        exit(8);
#endif
      }
    }
    Fhdf.close();
    */
    HDF5CloseFile(Fhdf);

    for(j=0, nbodies=0; j<nusetypes; j++) {
        k=usetypes[j];
        nbodies+=hdf_header_info.npartTotal[k];
        nbodies+=((long long)(hdf_header_info.npartTotalHW[k]) << 32);
    }
    return nbodies;

}
//@}



/// \name Get the number of hdf files per snapshot
//@{
inline Int_t HDF_get_nfiles(char *fname, int ptype)
{
    char buf[2000],buf1[2000],buf2[2000];
    sprintf(buf1,"%s.0.hdf5",fname);
    sprintf(buf2,"%s.hdf5",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    //H5File Fhdf;
    hid_t Fhdf;
    HDF_Group_Names hdf_gnames;
    //to store the groups, data sets and their associated data spaces
    //Attribute headerattribs;
    hid_t headerattribs;
    HDF_Header hdf_header_info;
    //buffers to load data
    int intbuff;
    long long longbuff;
    int ireaderror=0;
    Int_t nfiles = 0;
    //IntType inttype;

    //Try block to detect exceptions raised by any of the calls inside it
    //try
    {
        //turn off the auto-printing when failure occurs so that we can
        //handle the errors appropriately
        //Exception::dontPrint();

        //Open the specified file and the specified dataset in the file.
        //Fhdf.openFile(buf, H5F_ACC_RDONLY);
        Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
        //get header group
        hdf_header_info.num_files = read_attribute<int>(Fhdf, hdf_header_info.names[hdf_header_info.INumFiles]);
    }
    /*
    catch(GroupIException &error)
    {
        HDF5PrintError(error);
    }
    // catch failure caused by the H5File operations
    catch( FileIException &error )
    {
        HDF5PrintError(error);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException &error )
    {
        HDF5PrintError(error);
        ireaderror=1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException &error )
    {
        HDF5PrintError(error);
        ireaderror=1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException &error )
    {
        HDF5PrintError(error);
        ireaderror=1;
    }
    Fhdf.close();
    */
    HDF5CloseFile(Fhdf);

    return nfiles = hdf_header_info.num_files;

}
//@}

/// \name Wrappers to write attributes to HDF file
//@{
void WriteVELOCIraptorConfigToHDF(Options &opt, H5OutputFile &Fhdf);
///Write the simulation info (which could use input files to overwrite passed configuration options)
void WriteSimulationInfoToHDF(Options &opt, H5OutputFile &Fhdf);
///Write the unit info
void WriteUnitInfoToHDF(Options &opt, H5OutputFile &Fhdf);
//@}
#endif
