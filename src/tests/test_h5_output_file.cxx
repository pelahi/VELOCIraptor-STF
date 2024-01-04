#include <array>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef USEMPI
#include <mpi.h>
#endif // USEMPI

#include "allvars.h"
#include "hdfitems.h"
#include "logging.h"
#include "proto.h"

std::vector<double> generate_vector(std::size_t size)
{
    LOG(info) << "Creating vector of size " << size << " (" << size * sizeof(double) << " bytes)";
    std::vector<double> v(size);
    std::iota(v.begin(), v.end(), 0);
    return v;
}

template <typename T, typename U>
auto ceil_div(T a, U b) -> decltype(a + b)
{
    return (a + b - 1) / b;
}

int main(int argc, char *argv[])
{
#ifdef USEMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
#else
    int NProcs = 1;
    int ThisTask = 0;
#endif // USEMPI
    if (argc - 2 < NProcs) {
        std::cerr << "Usage: " << argv[0] << " <output-file> size1 size2 ... (one size per rank)\n";
        return 1;
    }

    std::string outfile {argv[1]};
    vr::init_logging(vr::LogLevel::trace);
    Options opts;
    opts.mpinprocswritesize = NProcs;
#ifdef USEMPI
    MPIInitWriteComm();
    MPIBuildWriteComm(opts);
#else
    int ThisWriteTask = 0;
#endif // USEMPI

    auto array_size = std::stoull(argv[ThisTask + 2]);
    H5OutputFile out;
    out.set_verbose(true);

    // Test serial writing of strings
    {
        out.create(outfile, H5F_ACC_TRUNC, 0, false);
        if (ThisWriteTask == 0) {
            out.write_dataset(opts, "string", 1, std::string("my example string"), false);
        }
        out.close();
    }

    {
        // Parallel writing of 1d array
        auto data = generate_vector(array_size);
        out.append(outfile);
        out.write_dataset(opts, "1d-doubles", data.size(), data.data());
        out.close();
    }

    {
        // 2 dimensions
        hsize_t dims2[2] = {ceil_div(array_size, 5), 5};
        auto data = generate_vector(dims2[0] * dims2[1]);
        out.append(outfile);
        out.write_dataset_nd(opts, "2d-doubles", 2, dims2, data.data());
        out.close();
    }
    {
        // 3 dimensions
        hsize_t dims3[3] = {ceil_div(array_size, 15), 3, 5};
        auto data = generate_vector(dims3[0] * dims3[1] * dims3[2]);
        out.append(outfile);
        out.write_dataset_nd(opts, "3d-doubles", 3, dims3, data.data());
        out.close();
    }

    // Attribute writing in serial
    {
        H5OutputFile out;
        out.set_verbose(true);
        out.append(outfile, H5F_ACC_RDWR, 0, false);
        if (ThisWriteTask == 0) {
            out.write_attribute("1d-doubles", "attribute_name", std::string("value"));
        }
        out.close();
    }

#ifdef USEMPI
    MPIFreeWriteComm();
    MPI_Finalize();
#endif // USEMPI
}