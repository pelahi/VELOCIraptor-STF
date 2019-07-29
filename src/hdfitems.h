/*! \file hdfitems.h
 *  \brief this file contains definitions and routines for reading HDF5 files
 *
 *   NOTE: the routines are based on reading Illustris HDF outputs
 */

#ifndef HDFITEMS_H
#define HDFITEMS_H


//#include "H5Cpp.h"
//using namespace H5;
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

#define HDFBHIMDOT 50
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
#define HDFNUMNAMETYPES  8
#define HDFILLUSTISNAMES 0
#define HDFGADGETXNAMES  1
#define HDFEAGLENAMES    2
#define HDFGIZMONAMES    3
#define HDFSIMBANAMES    4
#define HDFMUFASANAMES   5
#define HDFSWIFTEAGLENAMES    6
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
static inline hid_t hdf5_type(int dummy)                {return H5T_NATIVE_INT;}
static inline hid_t hdf5_type(long dummy)               {return H5T_NATIVE_LONG;}
static inline hid_t hdf5_type(long long dummy)          {return H5T_NATIVE_LLONG;}
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
	const hsize_t datarank, const hsize_t ndim, int nchunk, int noffset
)
{
	//setup hyperslab so that it is loaded into the buffer
	vector<hsize_t> start, count, stride, block, memdims;
	hid_t memspace;
	start.push_back(noffset);start.push_back(0);
	count.push_back(nchunk);count.push_back(ndim);
	stride.push_back(1);stride.push_back(1);
	block.push_back(1);block.push_back(1);
	memdims.push_back(nchunk*ndim);
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), stride.data(), count.data(), block.data());
	memspace = H5Screate_simple (1, memdims.data(), NULL);
	safe_hdf5<herr_t>(H5Dread, dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, buffer);

}

static inline void HDF5ReadHyperSlabInteger(long long *buffer,
	const hid_t &dataset, const hid_t &dataspace,
	const hsize_t datarank, const hsize_t ndim, int nchunk, int noffset
)
{
	//setup hyperslab so that it is loaded into the buffer
	vector<hsize_t> start, count, stride, block, memdims;
	hid_t memspace;
	start.push_back(noffset);start.push_back(0);
	count.push_back(nchunk);count.push_back(ndim);
	stride.push_back(1);stride.push_back(1);
	block.push_back(1);block.push_back(1);
	memdims.push_back(nchunk*ndim);
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start.data(), stride.data(), count.data(), block.data());
    memspace = H5Screate_simple (1, memdims.data(), NULL);
	safe_hdf5<herr_t>(H5Dread, dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT, buffer);
}

///\name HDF class to manage writing information
class H5OutputFile
{
	protected:

	hid_t file_id;

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
	}

	// Create a new file
	void create(std::string filename, unsigned int flag)
	{
		if(file_id >= 0)io_error("Attempted to create file when already open!");
		file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if(file_id < 0)io_error(string("Failed to create output file: ")+filename);
	}

    void append(std::string filename, unsigned int flag)
	{
		if(file_id >= 0)io_error("Attempted to create file when already open!");
		file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if(file_id < 0)io_error(string("Failed to create output file: ")+filename);
	}

	// Close the file
	void close()
	{
		if(file_id < 0)io_error("Attempted to close file which is not open!");
		H5Fclose(file_id);
		file_id = -1;
	}

  	// Destructor closes the file if it's open
	~H5OutputFile()
	{
	  if(file_id >= 0)
	    close();
	}

	/// Write a new 1D dataset. Data type of the new dataset is taken to be the type of
	/// the input data if not explicitly specified with the filetype_id parameter.
	template <typename T> void write_dataset(std::string name, hsize_t len, T *data,
	                                       hid_t memtype_id = -1, hid_t filetype_id=-1)
    {
		int rank = 1;
      	hsize_t dims[1] = {len};
		if (memtype_id == -1) {
    		memtype_id = hdf5_type(T{});
	 	}
      	write_dataset_nd(name, rank, dims, data, memtype_id, filetype_id);
    }
	void write_dataset(string name, hsize_t len, string data)
    {
		int rank = 1;
      	hsize_t dims[1] = {len};

		hid_t memtype_id, filetype_id, dspace_id, dset_id;
		herr_t status;
		memtype_id = H5Tcopy (H5T_C_S1);
		status = H5Tset_size (memtype_id, data.size());
		//status = H5Tset_size (memtype_id, H5T_VARIABLE);
		filetype_id = H5Tcopy (H5T_C_S1);
		//status = H5Tset_size (filetype_id, H5T_VARIABLE);
		status = H5Tset_size (filetype_id, data.size());

		// Create the dataspace
		dspace_id = H5Screate_simple(rank, dims, NULL);

		// Create the dataset
		dset_id = H5Dcreate(file_id, name.c_str(), filetype_id, dspace_id,
		                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		// Write the data
		if(H5Dwrite(dset_id, memtype_id, dspace_id, H5S_ALL, H5P_DEFAULT, data.c_str()) < 0)
		io_error(string("Failed to write dataset: ")+name);

		// Clean up (note that dtype_id is NOT a new object so don't need to close it)
		H5Sclose(dspace_id);
		H5Dclose(dset_id);
    }
	void write_dataset(string name, hsize_t len, void *data,
	                                       hid_t memtype_id=-1, hid_t filetype_id=-1)
    {
		int rank = 1;
      	hsize_t dims[1] = {len};
		if (memtype_id == -1) {
			throw std::runtime_error("Write data set called with void pointer but no type info passed.");
        }
      	write_dataset_nd(name, rank, dims, data, memtype_id, filetype_id);
    }


	/// Write a multidimensional dataset. Data type of the new dataset is taken to be the type of
	/// the input data if not explicitly specified with the filetype_id parameter.
	template <typename T> void write_dataset_nd(std::string name, int rank, hsize_t *dims, T *data,
	                                          hid_t memtype_id = -1, hid_t filetype_id=-1)
    {
		// Get HDF5 data type of the array in memory
		if (memtype_id == -1) {
			memtype_id = hdf5_type(T{});
		}

		// Determine type of the dataset to create
		if(filetype_id < 0)filetype_id = memtype_id;

		// Create the dataspace
		hid_t dspace_id = H5Screate_simple(rank, dims, NULL);

		// Only chunk non-zero size datasets
		int nonzero_size = 1;
		for(int i=0; i<rank; i+=1)
		if(dims[i]==0)nonzero_size = 0;

		// Only chunk datasets where we would have >1 chunk
		int large_dataset = 0;
		for(int i=0; i<rank; i+=1)
		if(dims[i] > HDFOUTPUTCHUNKSIZE)large_dataset = 1;

		// Dataset creation properties
		hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
		if(nonzero_size && large_dataset)
		{
			hsize_t *chunks = new hsize_t[rank];
			for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, dims[i]);
			H5Pset_layout(prop_id, H5D_CHUNKED);
			H5Pset_chunk(prop_id, rank, chunks);
			H5Pset_deflate(prop_id, HDFDEFLATE);
			delete chunks;
		}

		// Create the dataset
		hid_t dset_id = H5Dcreate(file_id, name.c_str(), filetype_id, dspace_id,
		                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(dset_id < 0)io_error(string("Failed to create dataset: ")+name);

		// Write the data
		if(H5Dwrite(dset_id, memtype_id, dspace_id, H5S_ALL, H5P_DEFAULT, data) < 0)
		io_error(string("Failed to write dataset: ")+name);

		// Clean up (note that dtype_id is NOT a new object so don't need to close it)
		H5Sclose(dspace_id);
		H5Dclose(dset_id);
		H5Pclose(prop_id);

    }
	void write_dataset_nd(std::string name, int rank, hsize_t *dims, void *data,
	                                          hid_t memtype_id = -1, hid_t filetype_id=-1)
    {
		// Get HDF5 data type of the array in memory
		if (memtype_id == -1) {
			throw std::runtime_error("Write data set called with void pointer but no type info passed.");
		}

		// Determine type of the dataset to create
		if(filetype_id < 0)filetype_id = memtype_id;

		// Create the dataspace
		hid_t dspace_id = H5Screate_simple(rank, dims, NULL);

		// Only chunk non-zero size datasets
		int nonzero_size = 1;
		for(int i=0; i<rank; i+=1)
		if(dims[i]==0)nonzero_size = 0;

		// Only chunk datasets where we would have >1 chunk
		int large_dataset = 0;
		for(int i=0; i<rank; i+=1)
		if(dims[i] > HDFOUTPUTCHUNKSIZE)large_dataset = 1;

		// Dataset creation properties
		hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
		if(nonzero_size && large_dataset)
		{
			hsize_t *chunks = new hsize_t[rank];
			for(auto i=0; i<rank; i++) chunks[i] = min((hsize_t) HDFOUTPUTCHUNKSIZE, dims[i]);
			H5Pset_layout(prop_id, H5D_CHUNKED);
			H5Pset_chunk(prop_id, rank, chunks);
			H5Pset_deflate(prop_id, HDFDEFLATE);
			delete chunks;
		}

		// Create the dataset
		hid_t dset_id = H5Dcreate(file_id, name.c_str(), filetype_id, dspace_id,
		                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(dset_id < 0)io_error(string("Failed to create dataset: ")+name);

		// Write the data
		if(H5Dwrite(dset_id, memtype_id, dspace_id, H5S_ALL, H5P_DEFAULT, data) < 0)
		io_error(string("Failed to write dataset: ")+name);

		// Clean up (note that dtype_id is NOT a new object so don't need to close it)
		H5Sclose(dspace_id);
		H5Dclose(dset_id);
		H5Pclose(prop_id);

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
		hid_t attr_id = H5Acreate(file_id, name.c_str(), dtype_id, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
		if(attr_id < 0)io_error(string("Unable to create attribute ")+name+string(" on object ")+parent);

		// Write the attribute
		if(H5Awrite(attr_id, dtype_id, &data) < 0)
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
            names[itemp++]=string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");
            names[itemp++]=string("ParticleIDs");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Mass");
            else names[itemp++]=string("Masses");
            names[itemp++]=string("Density");
            names[itemp++]=string("InternalEnergy");
            names[itemp++]=string("StarFormationRate");
            //always place the metacallity at position 7 in naming array
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
            else if (hdfnametype==HDFEAGLENAMES || hdfnametype==HDFSWIFTEAGLENAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=string("Metallicity");
            }
        }
        //dark matter
        if (ptype==HDFDMTYPE) {
            names[itemp++]=string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");
            names[itemp++]=string("ParticleIDs");
            if (hdfnametype==HDFSWIFTEAGLENAMES) {
                names[itemp++]=string("Masses");
            }
            if (hdfnametype==HDFSIMBANAMES||hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=string("Masses");
                names[itemp++]=string("Potential");
            }
            if (hdfnametype==HDFILLUSTISNAMES) {
                names[itemp++]=string("Potential");
                names[itemp++]=string("SubfindDensity");
                names[itemp++]=string("SubfindHsml");
                names[itemp++]=string("SubfindVelDisp");
            }
        }
        //also dark matter particles
        if (ptype==HDFDM1TYPE ||ptype==HDFDM2TYPE) {
            names[itemp++]=string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");
            names[itemp++]=string("ParticleIDs");
            names[itemp++]=string("Masses");
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
            names[itemp++]=string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");
            names[itemp++]=string("ParticleIDs");
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
        }
        if (ptype==HDFBHTYPE) {
            names[itemp++]=string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Velocity");
            else names[itemp++]=string("Velocities");
            names[itemp++]=string("ParticleIDs");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=string("Mass");
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

        if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) {

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

}//@}



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

#endif
