#ifndef H5_UTILS_H
#define H5_UTILS_H

#ifdef USEHDF

#include <string>

#include <hdf5.h>

namespace vr
{

/**
 * The function invoked when an error happens when invoking an HDF5 function.
 * @param fname The name of the faulty HDF5 function
 */
void _hdf5_error_handler(const std::string &fname);

/**
 * Invokes an HDF5 function with the given arguments, and automatically invoke
 * our own error handler in case of an error. Our handler aborts the program.
 *
 * This is meant to be invoked via the safe_hdf5 macro.
 *
 * @tparam F The function type
 * @tparam Args The argument types
 * @param function The HDF5 function to invoke
 * @param fname The HDF5 funciton name
 * @param args The arguments to invoke the HDF5 function with
 * @return The result of the HDF5 function invocation
 */
template <typename F, typename ... Args>
auto _safe_hdf5(F function, const std::string &fname, Args ... args) -> decltype(function(std::forward<Args>(args)...))
{
    using ReturnT = decltype(function(std::forward<Args>(args)...));
    ReturnT ret = function(std::forward<Args>(args)...);
    if (ret < 0) {
        _hdf5_error_handler(fname);
    }
    return ret;
}

/// Wrapper around _safe_hdf5, this is what users should invoke
#define safe_hdf5(F, ...) ::vr::_safe_hdf5(F, # F, __VA_ARGS__)

/**
 * Get a string with information about the given dataset.
 * @param dataspace_id The dataset ID
 * @return A string with information about the given dataset.
 */
std::string dataspace_information(hid_t dataspace_id);

} // namespace vr

#endif // USEHDF

#endif /* H5_UTILS_H */