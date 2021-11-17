#ifdef USEHDF

#include <sstream>
#include <vector>

#include "h5_utils.h"
#include "ioutils.h"
#include "logging.h"

#ifdef USEMPI
#include <mpi.h>
#endif // USEMPI

namespace vr
{

void _hdf5_error_handler(const std::string &fname)
{
    LOG(error) << "Error in HDF routine " << fname << ", aborting";
#ifdef USEMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    abort();
#endif
}

std::string dataspace_information(hid_t dataspace_id)
{
    std::ostringstream os;
    int ndims = H5Sget_simple_extent_ndims(dataspace_id);
    os << "Dataspace " << dataspace_id << " has " << ndims << " dimensions. ";

    std::vector<hsize_t> dims(ndims);
    H5Sget_simple_extent_dims(dataspace_id, dims.data(), nullptr);
    os << "Size: " << vr::printable_range(dims) << ". ";

    auto selection_type = H5Sget_select_type(dataspace_id);
    if (selection_type == H5S_SEL_HYPERSLABS) {

#if H5_VERSION_GE(1, 10, 0)
        if (safe_hdf5(H5Sis_regular_hyperslab, dataspace_id)) {
            std::vector<hsize_t> start(ndims);
            std::vector<hsize_t> stride(ndims);
            std::vector<hsize_t> count(ndims);
            std::vector<hsize_t> block(ndims);
            safe_hdf5(H5Sget_regular_hyperslab, dataspace_id, start.data(), stride.data(), count.data(), block.data());
            os << "Hyperslab: start=" << vr::printable_range(start)
               << ", stride=" << vr::printable_range(stride)
               << ", count=" << vr::printable_range(count)
               << ", block=" << vr::printable_range(block);
        }
        else
#endif
        {
            hssize_t nhyperblocks = H5Sget_select_hyper_nblocks(dataspace_id);
            std::vector<hsize_t> hyperblocks(nhyperblocks * ndims * 2);
            H5Sget_select_hyper_blocklist(dataspace_id, 0, nhyperblocks, hyperblocks.data());

            // Create "slabs" from start-end "blocks"
            std::vector<std::pair<std::vector<hsize_t>, std::vector<hsize_t>>> slabs(nhyperblocks);
            auto coord_start = hyperblocks.begin();
            for (std::size_t block = 0; block != nhyperblocks; block++) {
                auto coord_end = coord_start + ndims;
                slabs[block].first = std::vector<hsize_t>(coord_start, coord_end);
                coord_start += ndims * 2;
                slabs[block].second = std::vector<hsize_t>(coord_end, coord_start);
            }
            os << "Hyperslabs: ";
            for (const auto &slab : slabs) {
                os << vr::printable_range(slab.first);
                os << " -> ";
                os << vr::printable_range(slab.second);
                os << " ";
            }
        }

    }
    else if (selection_type == H5S_SEL_NONE) {
        os << "Dataspace fully deselected";
    }
    else if (selection_type == H5S_SEL_ALL) {
        os << "Dataspace fully selected";
    }
    else if (selection_type == H5S_SEL_POINTS) {
        os << "Dataspace selected by points, not shown here";
    }
    return os.str();
};

}  // namespace vr

#endif // USEHDF