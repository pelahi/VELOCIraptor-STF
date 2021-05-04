#include <sstream>

#ifdef USEOPENMP
#include <omp.h>
#endif // USEOPENMP

#include "exceptions.h"

#include "allvars.h"

namespace vr
{

static std::string build_error_message(const Particle &part, const std::string &function)
{
    std::ostringstream os;
    os << "Particle density not positive, cannot continue.\n\n"
       << "Particle information: "
       << "id=" << part.GetID() << ", pid=" << part.GetPID()
       << ", type=" << part.GetType()
       << ", pos=(" << part.X() << ", " << part.Y() << ", " << part.Z()
       << "), vel=(" << part.Vx() << ", " << part.Vy() << ", " << part.Vz()
       << "), mass=" << part.GetMass()
       << ", density=" << part.GetDensity() << '\n';
    os << "Execution context: MPI enabled=";
#ifdef USEMPI
    os << "yes, rank=" << ThisTask << '/' << NProcs;
#else
    os << "no";
#endif // USEMPI
    os << ", OpenMP enabled=";
#ifdef USEOPENMP
    os << "yes, thread=" << omp_get_thread_num() << '/' << omp_get_num_threads();
#else
    os << "no";
#endif // USEOPENMP
    return os.str();
}

non_positive_density::non_positive_density(const Particle &part, const std::string &function)
  : exception(build_error_message(part, function))
{
}

}  // namespace vr