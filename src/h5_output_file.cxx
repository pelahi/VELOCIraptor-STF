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