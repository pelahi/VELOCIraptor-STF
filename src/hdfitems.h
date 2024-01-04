/*! \file hdfitems.h
 *  \brief this file contains definitions and routines for reading HDF5 files
 *
 *   NOTE: the routines are based on reading Illustris HDF outputs
 */

#ifndef HDFITEMS_H
#define HDFITEMS_H

#include <hdf5.h>
#include <string>

#include "h5_utils.h"
#include "ioutils.h"
#include "io.h"
#include "logging.h"


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
#define HDFHEADNINFO 20
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
#define HDFGASTEMP 99
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
#define HDFSWIFTFLAMINGONAMES    9
//@}

///size of chunks in hdf files for Compression
#define HDFOUTPUTCHUNKSIZE 8192
#define HDFDEFLATE    6


#if H5_VERSION_GE(1,10,1)
#define HDF5_FILE_GROUP_COMMON_BASE H5::Group
#else
#define HDF5_FILE_GROUP_COMMON_BASE H5::CommonFG
#endif

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

static inline int whatisopen(hid_t fid) {
    ssize_t cnt;
    int howmany;
    int i;
    H5I_type_t ot;
    hid_t anobj;
    hid_t *objs;
    char name[1024];
    herr_t status;
    cnt = safe_hdf5(H5Fget_obj_count, fid, H5F_OBJ_ALL);
    LOG(info) << cnt << " object(s) open";
    if (cnt <= 0) return cnt;
    objs = new hid_t[cnt];
    howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, cnt, objs);
    LOG(info) << "open objects:";
    for (i = 0; i < howmany; i++ ) {
            anobj = objs[i];
            ot = H5Iget_type(anobj);
            status = H5Iget_name(anobj, name, 1024);
            LOG(info) << i <<" type "<<ot<<" name "<<name;
    }
    delete[] objs;
    return howmany;
}

//template <typename AttributeHolder>
//static inline H5::Attribute get_attribute(const AttributeHolder &l, const std::string attr_name)
static inline void get_attribute(vector<hid_t> &ids, const std::string attr_name)
{
    //can use H5Aexists as it is the C interface but how to access it?
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
        hid_t lapl_id = H5P_DEFAULT;
#if H5_VERSION_GE(1,12,0)
        unsigned int fields = H5O_INFO_ALL;
        H5Oget_info_by_name(ids.back(), parts[0].c_str(), &object_info, fields, lapl_id);
#else
        H5Oget_info_by_name(ids.back(), parts[0].c_str(), &object_info, lapl_id);
#endif
        if (object_info.type == H5O_TYPE_GROUP) {
            newid = H5Gopen(ids.back(),parts[0].c_str(),H5P_DEFAULT);
        }
        else if (object_info.type == H5O_TYPE_DATASET) {
            newid = H5Dopen(ids.back(),parts[0].c_str(),H5P_DEFAULT);
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
        hid_t lapl_id = H5P_DEFAULT;
#if H5_VERSION_GE(1,12,0)
        H5Oget_info(id, &object_info, lapl_id);
#else
        H5Oget_info(id, &object_info);
#endif
        auto ot = H5Iget_type(id);
        if (ot == H5I_GROUP) {
            safe_hdf5(H5Gclose, id);
        }
        else if (ot == H5I_DATASET) {
            safe_hdf5(H5Dclose, id);
        }
        else if (ot == H5I_DATASPACE) {
            safe_hdf5(H5Sclose, id);
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

template<typename T> T read_attribute(const hid_t &file_id, const std::string &name) {
    std::string attr_name;
    T val;
    hid_t type;
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
template<typename T> vector<T> read_attribute_v(const hid_t &file_id, const std::string &name) {
    std::string attr_name;
    vector<T> val;
    hid_t type;
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

static std::string get_hdf5_name(const hid_t object_id)
{
    char name[128];
    if (H5Iget_name(object_id, name, 128) == 0) {
        throw std::runtime_error("No name for object " + std::to_string(object_id));
    }
    return name;
}

template<typename T> T read_attribute(const std::string &filename, const std::string &name) {
    LOG_RANK0(debug) << "Reading attribute " << name;
    safe_hdf5(H5Fopen, filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    T attr = read_attribute<T>(file_id, name);
    safe_hdf5(H5Fclose,file_id);
    return attr;
}

static inline hid_t HDF5OpenFile(string name, unsigned int flags){
    LOG_RANK0(info) << "Opening " << name;
    return H5Fopen(name.c_str(),flags, H5P_DEFAULT);
}
static inline hid_t HDF5OpenGroup(const hid_t &file, string name){
    LOG_RANK0(debug) << "Opening Group " << name << '/' << name;
    hid_t idval = safe_hdf5(H5Gopen, file, name.c_str(), H5P_DEFAULT);
    return idval;
}
static inline hid_t HDF5OpenDataSet(const hid_t &id, string name){
    LOG_RANK0(debug) << "Opening Dataset " << get_hdf5_name(id) << '/' << name;
    hid_t idval = safe_hdf5(H5Dopen, id, name.c_str(), H5P_DEFAULT);
    return idval;
}
static inline hid_t HDF5OpenDataSpace(const hid_t &id){
    hid_t idval=H5Dget_space(id);
    return idval;
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
static inline void HDF5CloseAll(hid_t &fid){
    ssize_t cnt;
    std::vector<hid_t> objs;
    cnt = safe_hdf5(H5Fget_obj_count, fid, H5F_OBJ_ALL);
    LOG(info) <<"values are "<<cnt;
    if (cnt <= 0) return;
    objs.resize(cnt);
    safe_hdf5(H5Fget_obj_ids, fid, H5F_OBJ_ALL, cnt, objs.data());
    close_hdf_ids(objs);
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
    safe_hdf5(H5Dread, dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, plist_id, buffer);

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
    safe_hdf5(H5Dread, dataset, H5T_NATIVE_LONG, memspace, dataspace, plist_id, buffer);
}

///\name HDF class to manage writing information
class H5OutputFile
{
public:
    constexpr static int ALL_RANKS = -1;

private:
    bool verbose = false;
    hid_t file_id = -1;
    int writing_rank = ALL_RANKS;

    // Called if a HDF5 call fails (might need to MPI_Abort)
    void io_error(std::string message) {
        std::cerr << message << std::endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
        abort();
    }

    void truncate(const std::string &filename, hid_t access_plist);
    void open(const std::string &filename, hid_t access_plist);

    template <typename CreationFunction>
    void create(CreationFunction file_creator, hid_t flag, int rank);

    void write_attribute(const std::string &parent, const std::string &name, hid_t dtype_id, const void *data);

public:

    void set_verbose(bool verbose)
    {
        this->verbose = verbose;
    }

    // Create a new file
    void create(std::string filename, hid_t flag = H5F_ACC_TRUNC,
        int taskID = ALL_RANKS, bool iparallelopen = true);

    void append(std::string filename, hid_t flag = H5F_ACC_RDWR,
        int taskID = ALL_RANKS, bool iparallelopen = true);

    // Close the file
    void close();

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
    template <typename T> void write_dataset(Options opt, std::string name, hsize_t len, T *data,
       hid_t memtype_id = -1, hid_t filetype_id=-1, bool flag_parallel = true)
    {
        assert(memtype_id == -1);
        assert(filetype_id == -1);
        memtype_id = hdf5_type(T{});
        write_dataset(opt, name, len, static_cast<void *>(data),
                      memtype_id, filetype_id, flag_parallel);
    }

    void write_dataset(Options opt, string name, hsize_t len, void *data,
       hid_t memtype_id, hid_t filetype_id = -1, bool flag_parallel = true);

    void write_dataset(Options opt, string name, hsize_t len, string data, bool flag_parallel = true);

    /// Write a multidimensional dataset. Data type of the new dataset is taken to be the type of
    /// the input data if not explicitly specified with the filetype_id parameter. This is just
    /// a wrapper which uses the type of the supplied array to set the memtype_id parameter
    /// if it's not set explicitly.
    template <typename T> void write_dataset_nd(Options opt, std::string name, int rank, hsize_t *dims, T *data,
        hid_t memtype_id = -1, hid_t filetype_id = -1,
        bool flag_parallel = true)
    {
        if (memtype_id == -1) {
            memtype_id = hdf5_type(T{});
        }
        write_dataset_nd(opt, name, rank, dims, (void *) data,
                         memtype_id, filetype_id, flag_parallel);
    }

    void write_dataset_nd(Options opt, std::string name, int rank, hsize_t *dims, void *data,
        hid_t memtype_id = -1, hid_t filetype_id = -1, bool flag_parallel = true);

    /// write an attribute
    template <typename T> void write_attribute(std::string parent, std::string name, T data)
    {
        write_attribute(parent, name, hdf5_type(data), &data);
    }

    void write_attribute(string parent, string name, string data);
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
	  case HDFSWIFTFLAMINGONAMES:
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
    double Omegab, Omegacdm, Omegar, Omegade, Omegak;
    double wde, wde0, wdea;

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
    const static int IOmegab   =12;
    const static int IOmegar   =13;
    const static int IOmegak   =14;
    const static int Iwde      =15;
    const static int Iwde0     =16;
    const static int Iwdea     =17;

    ///constructor
    HDF_Header(int hdfnametype=HDFEAGLENAMES) {
        int itemp=0;
        switch (hdfnametype) {
          case HDFSWIFTEAGLENAMES:
    	  case HDFSWIFTFLAMINGONAMES:
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
            names[itemp++]=string("Cosmology/Omega_b");
            names[itemp++]=string("Cosmology/Omega_r");
            names[itemp++]=string("Cosmology/Omega_k");
            names[itemp++]=string("Cosmology/w");
            names[itemp++]=string("Cosmology/w_0");
            names[itemp++]=string("Cosmology/w_a");
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

    //the HDF naming convention for the data blocks. By default assumes ILLUSTRIS naming convention
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
	    else if(hdfnametype==HDFSWIFTFLAMINGONAMES) names[itemp++]=string("Densities");
            else names[itemp++]=string("Density");

            // Internal energies
            if(hdfnametype==HDFSWIFTEAGLENAMES) names[itemp++]=string("InternalEnergies");
            else names[itemp++]=string("InternalEnergy");

            // SFR
            if(hdfnametype==HDFSWIFTEAGLENAMES) names[itemp++]=string("StarFormationRates");
	    else if(hdfnametype==HDFSWIFTFLAMINGONAMES) names[itemp++]=string("StarFormationRates");
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
            else if(hdfnametype==HDFSWIFTEAGLENAMES || hdfnametype==HDFSWIFTFLAMINGONAMES) {
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
                hdfnametype==HDFMUFASANAMES || hdfnametype==HDFOLDSWIFTEAGLENAMES ||
		hdfnametype==HDFSWIFTFLAMINGONAMES) {
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
            else if (hdfnametype==HDFSWIFTEAGLENAMES ||
		     hdfnametype==HDFSWIFTFLAMINGONAMES) {
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
            else if(hdfnametype==HDFSWIFTEAGLENAMES || hdfnametype==HDFSWIFTFLAMINGONAMES)
	      names[itemp++]=string("DynamicalMasses");
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
            else if (hdfnametype==HDFSWIFTEAGLENAMES ||
		     hdfnametype==HDFSWIFTFLAMINGONAMES) {
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
	if (ptype==HDFGASTYPE) {
            // Temperature
            propindex[HDFGASTEMP] = itemp;
            if(hdfnametype==HDFSWIFTEAGLENAMES || hdfnametype==HDFSWIFTFLAMINGONAMES)
	      names[itemp++]=string("Temperatures");
            else names[itemp++]=string("Temperature");
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
            if (opt.ihdfnameconvention!=HDFSWIFTEAGLENAMES && opt.ihdfnameconvention!=HDFSWIFTFLAMINGONAMES)
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
            if (opt.ihdfnameconvention!=HDFSWIFTEAGLENAMES && opt.ihdfnameconvention!=HDFSWIFTFLAMINGONAMES)
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

static std::string find_hdf5_file(const char *prefix)
{
    for (auto *suffix : {".0.hdf5", ".hdf5"}) {
        std::string filename(prefix);
        filename += suffix;
        LOG(debug) << "Looking for " << filename;
        if (FileExists(filename.c_str())) {
            return filename;
        }
    }
    LOG(error) << "Can't find HDF5 file with prefix " << prefix;
    exit(9);
}

/// \name Get the number of particles in the hdf files
//@{
inline Int_t HDF_get_nbodies(char *fname, int ptype, Options &opt)
{
    auto buf = find_hdf5_file(fname);

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
    vector<unsigned int> vuintbuff;
    int j,k;
    Int_t nbodies=0;

    //to determine types
    //IntType inttype;
    //StrType stringtype;
    int nusetypes,usetypes[NHDFTYPE],nbusetypes;
    HDFSetUsedParticleTypes(opt,nusetypes,nbusetypes,usetypes);

    Fhdf = H5Fopen(buf.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    LOG(info) << "Loading HDF header info in header group: " << hdf_gnames.Header_name;

    if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFSWIFTFLAMINGONAMES) {

        // Check if it is a SWIFT snapshot.
        dataname = string("Header/Code");
        stringbuff = read_attribute<string>(Fhdf, dataname);

        // Read SWIFT parameters
        if(!swift_str.compare(stringbuff)) {
            hdf_header_info.iscosmological=read_attribute<int>(Fhdf, hdf_header_info.names[hdf_header_info.IIsCosmological]);
            if (!hdf_header_info.iscosmological && opt.icosmologicalin) {
                LOG(error)<<"Error: cosmology is turned on in the config file but the snaphot provided is a non-cosmological run.";
#ifdef USEMPI
                MPI_Abort(MPI_COMM_WORLD, 8);
#else
                exit(0);
#endif
            }
            else if (hdf_header_info.iscosmological && !opt.icosmologicalin) {
                LOG(error)<<"Error: cosmology is turned off in the config file but the snaphot provided is a cosmological run.";
#ifdef USEMPI
                MPI_Abort(MPI_COMM_WORLD, 8);
#else
                exit(0);
#endif
            }
            //swift does not have little h's in input so don't convert
            opt.inputcontainslittleh = false;

        }
        // If the code is not SWIFT
        else {
            LOG(error)<<"SWIFT EAGLE HDF5 naming convention chosen in config file but the snapshot was not produced by SWIFT. The string read was: "<<stringbuff;
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
            LOG(warning)<<"Configured to load star particles but none present";
            opt.iusestarparticles=0;
        }
        if (opt.iusesinkparticles && hdf_header_info.npartTotalHW[HDFBHTYPE] == 0 && hdf_header_info.npartTotal[HDFBHTYPE] == 0)
        {
            LOG(warning)<<"Configured to load black hole particles but none present";
            opt.iusesinkparticles=0;
        }
        if (opt.iusewindparticles && hdf_header_info.npartTotalHW[HDFWINDTYPE] == 0 && hdf_header_info.npartTotal[HDFWINDTYPE] == 0)
        {
            LOG(warning)<<"Configured to load wind particles but none present";
            opt.iusewindparticles=0;
        }
        if (opt.iusetracerparticles && hdf_header_info.npartTotalHW[HDFTRACERTYPE] == 0 && hdf_header_info.npartTotal[HDFTRACERTYPE] == 0)
        {
            LOG(warning)<<"Configured to load tracer particles but none present";
            opt.iusetracerparticles=0;
        }
#ifdef HIGHRES
        if (opt.iuseextradarkparticles && hdf_header_info.npartTotalHW[HDFDM1TYPE] == 0 
            && hdf_header_info.npartTotal[HDFDM1TYPE] == 0
            && hdf_header_info.npartTotalHW[HDFDM2TYPE] == 0 
            && hdf_header_info.npartTotal[HDFDM2TYPE] == 0)
        {
            LOG(warning)<<"Configured to load extra low res dark matter particles but none present";
            opt.iuseextradarkparticles=0;
        }
#endif
    }
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
    auto buf = find_hdf5_file(fname);

    //H5File Fhdf;
    hid_t Fhdf;
    HDF_Group_Names hdf_gnames;
    //to store the groups, data sets and their associated data spaces
    //Attribute headerattribs;
    hid_t headerattribs;
    HDF_Header hdf_header_info;
    Int_t nfiles = 0;
    //IntType inttype;


    //Open the specified file and the specified dataset in the file.
    //Fhdf.openFile(buf, H5F_ACC_RDONLY);
    Fhdf = H5Fopen(buf.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    //get header group
    hdf_header_info.num_files = read_attribute<int>(Fhdf, hdf_header_info.names[hdf_header_info.INumFiles]);
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
