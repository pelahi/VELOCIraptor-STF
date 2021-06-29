/*! \file h5_output_file.cxx
 */

#include "hdfitems.h"
#include "io.h"

void H5OutputFile::truncate(const std::string &filename, hid_t access_plist)
{
    file_id = safe_hdf5(H5Fcreate, filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access_plist);
}

void H5OutputFile::open(const std::string &filename, hid_t access_plist)
{
    file_id = safe_hdf5(H5Fopen, filename.c_str(), H5F_ACC_RDWR, access_plist);
}

template <typename CreationFunction>
void H5OutputFile::create(CreationFunction file_creator, hid_t flag, int rank)
{
    writing_rank = rank;
    if (file_id >= 0) {
        io_error("Attempted to create file when already open!");
    }
#ifdef USEPARALLELHDF
    MPI_Comm comm = mpi_comm_write;
    if (writing_rank == ALL_RANKS) {
        auto parallel_access_id = safe_hdf5(H5Pcreate, H5P_FILE_ACCESS);
        safe_hdf5(H5Pset_fapl_mpio, parallel_access_id, comm, MPI_INFO_NULL);
        file_creator(parallel_access_id);
        safe_hdf5(H5Pclose, parallel_access_id);
    }
    else {
        if (writing_rank >= NProcsWrite) {
            io_error(string("MPI rank asked to create file out of range. Rank is ") + to_string(writing_rank));
        }
        else if (ThisWriteTask == writing_rank) {
            file_creator(H5P_DEFAULT);
        }
        MPI_Barrier(comm); // TODO: check if really needed
    }
#else
    file_creator(H5P_DEFAULT);
#endif
}

void H5OutputFile::create(std::string filename, hid_t flag, int rank, bool iparallelopen)
{
    assert(iparallelopen == (rank == ALL_RANKS));
    assert(flag == H5F_ACC_TRUNC);
    create(
        [&](hid_t access_plist) {
            truncate(filename, access_plist);
        },
        flag, rank);
}

void H5OutputFile::append(std::string filename, hid_t flag, int rank, bool iparallelopen)
{
    assert(iparallelopen == (rank == ALL_RANKS));
    assert(flag == H5F_ACC_RDWR);
    create(
        [&](hid_t access_plist) {
            open(filename, access_plist);
        },
        flag, rank);
}

void H5OutputFile::close()
{
#ifdef USEPARALLELHDF
    // we didn't open anything for writing
    if (writing_rank != ALL_RANKS && writing_rank != ThisWriteTask) {
        return;
    }
#endif // USEPARALLELHDF
    if(file_id < 0) {
        io_error("Attempted to close file which is not open!");
    }
    H5Fclose(file_id);
    file_id = -1;
}


void H5OutputFile::write_dataset(Options opt, string name, hsize_t len, void *data,
   hid_t memtype_id, hid_t filetype_id, bool flag_parallel)
{
    int rank = 1;
    hsize_t dims[1] = {len};
    write_dataset_nd(opt, name, rank, dims, data, memtype_id, filetype_id, flag_parallel);
}

void H5OutputFile::write_dataset(Options opt, string name, hsize_t len, string data,
    bool flag_parallel)
{
#ifdef USEPARALLELHDF
    assert(!flag_parallel);
#endif
    int rank = 1;
    hsize_t dims[1] = {len};
    auto memtype_id = safe_hdf5(H5Tcopy, H5T_C_S1);
    safe_hdf5(H5Tset_size, memtype_id, data.size());
    auto filetype_id = safe_hdf5(H5Tcopy, H5T_C_S1);
    safe_hdf5(H5Tset_size, filetype_id, data.size());
    write_dataset_nd(opt, name, rank, dims, data.data(), memtype_id, filetype_id, false);
    safe_hdf5(H5Tclose, filetype_id);
    safe_hdf5(H5Tclose, memtype_id);
}


static hid_t get_dataset_creation_property(int rank, hsize_t *dims)
{
#ifdef USEHDFCOMPRESSION
    // Compress data only if all dimensions are > 0
    // and at least one dimension must be "large enough"
    auto positive_dims = std::all_of(dims, dims + rank, [](hsize_t dim) { return dim > 0; });
    auto large_dataset = std::any_of(dims, dims + rank, [](hsize_t dim) { return dim > HDFOUTPUTCHUNKSIZE; });

    if (positive_dims && large_dataset) {
        vector<hsize_t> chunks(rank);
        std::transform(dims, dims + rank, chunks.begin(),
            [](hsize_t dim) { return std::min(dim, static_cast<hsize_t>(HDFOUTPUTCHUNKSIZE)); }
        );
        auto prop_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(prop_id, H5D_CHUNKED);
        H5Pset_chunk(prop_id, rank, chunks.data());
        H5Pset_deflate(prop_id, HDFDEFLATE);
        return prop_id;
    }
#endif
    return H5P_DEFAULT;
}

void H5OutputFile::write_dataset_nd(Options opt, std::string name,
    int ndims, hsize_t *dims, void *data,
    hid_t memtype_id, hid_t filetype_id,
    bool flag_parallel)
{
    bool write_in_parallel = flag_parallel && opt.mpinprocswritesize > 1;
    assert(ndims > 0);

    if (verbose && LOG_ENABLED(debug)) {
        std::ostringstream os;
        os << "Writing dataset " << name << " with dimensions ";
        os << vr::printable_range(dims, ndims);
        os << " in " << (write_in_parallel ? "parallel" : "serial");
        LOG(debug) << os.str();
    }
    // Get HDF5 data type of the array in memory
    if (memtype_id == -1) {
        throw std::runtime_error("Write data set called with void pointer but no type info passed.");
    }
    // Determine type of the dataset to create
    if (filetype_id < 0) {
        filetype_id = memtype_id;
    }

    // Get full extent of the dataset and calculate local offsets, if required.
    // Only the first dimension is written in parallel, it's assumed
    // all other dimensions are the same in all MPI ranks
#ifdef USEPARALLELHDF
    std::vector<hsize_t> extended_dims(dims, dims + ndims);
    std::vector<hsize_t> dims_offset(ndims);
    if (write_in_parallel) {

        int comm_size;
        int comm_rank;
        MPI_Comm_size(mpi_comm_write, &comm_size);
        MPI_Comm_rank(mpi_comm_write, &comm_rank);

        std::vector<hsize_t> all_dims(ndims * comm_size);
        static_assert(sizeof(hsize_t) == sizeof(unsigned long long), "hsize_t not same size of unsigned long long");

        MPI_Allgather(
            dims, ndims, MPI_UNSIGNED_LONG_LONG,
            all_dims.data(), ndims, MPI_UNSIGNED_LONG_LONG,
            mpi_comm_write);

        // All dimensions other than the first one should be the same in all ranks
#ifndef NDEBUG
        for (auto rank = 0; rank != comm_size; rank++) {
            auto rank_dims_start = rank * ndims;
            for (auto dim = 1; dim != ndims; dim++) {
                assert(all_dims[dim] == all_dims[dim + rank_dims_start]);
            }
        }
#endif

        // Only first dimension is parallel and gets adjusted
        extended_dims[0] = 0;
        for (auto rank = 0; rank != comm_size; rank++) {
            auto first_dim_size = all_dims[rank * ndims];
            extended_dims[0] += first_dim_size;
            if (rank < comm_rank) {
                dims_offset[0] += first_dim_size;
            }
        }
    }
#endif

    // Create dataspaces. When writing in parallel, the file dataspace spans
    // the full extent of the (distributed) data, so it's different
    hid_t memspace_id = H5Screate_simple(ndims, dims, NULL);
    hid_t dspace_id = memspace_id;
#ifdef USEPARALLELHDF
    if (write_in_parallel) {
        dspace_id = H5Screate_simple(ndims, extended_dims.data(), NULL);
    }
#endif // USEPARALLELHDF

    // Create the dataset
    hid_t prop_id;
#ifdef USEPARALLELHDF
    if (write_in_parallel) {
        prop_id = get_dataset_creation_property(ndims, extended_dims.data());
    }
    else
#endif
    {
        prop_id = get_dataset_creation_property(ndims, dims);
    }
    auto dset_id = safe_hdf5(H5Dcreate, file_id, name.c_str(), filetype_id, dspace_id,
        H5P_DEFAULT, prop_id, H5P_DEFAULT);
    H5Pclose(prop_id);

    // Do the actual data writing
    prop_id = H5P_DEFAULT;
    auto write_data = [&]() {
        if (verbose) {
            LOG(debug) << vr::dataspace_information(dspace_id);
            LOG(debug) << vr::dataspace_information(memspace_id);
        }
        safe_hdf5(H5Dwrite, dset_id, memtype_id, memspace_id, dspace_id, prop_id, data);
    };
#ifdef USEPARALLELHDF
    if (write_in_parallel) {
        // set up the collective transfer properties list
        prop_id = safe_hdf5(H5Pcreate, H5P_DATASET_XFER);
        safe_hdf5(H5Pset_dxpl_mpio, prop_id, H5FD_MPIO_COLLECTIVE);
        safe_hdf5(H5Sselect_hyperslab, dspace_id, H5S_SELECT_SET, dims_offset.data(), nullptr, dims, nullptr);
        if (dims[0] == 0) {
            safe_hdf5(H5Sselect_none, dspace_id);
            safe_hdf5(H5Sselect_none, memspace_id);
        }
        if (extended_dims[0] > 0) {
            write_data();
        }
    }
    else
#endif // USEPARALLELHDF
    {
        if (dims[0] > 0) {
            write_data();
        }
    }

    // Clean up (note that dtype_id is NOT a new object so don't need to close it)
    safe_hdf5(H5Pclose, prop_id);
    if (memspace_id != dspace_id) {
        safe_hdf5(H5Sclose, memspace_id);
    }
    safe_hdf5(H5Sclose, dspace_id);
    safe_hdf5(H5Dclose, dset_id);
}