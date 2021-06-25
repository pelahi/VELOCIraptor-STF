#ifdef USEHDF

#include "h5_utils.h"
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

}  // namespace vr

#endif // USEHDF