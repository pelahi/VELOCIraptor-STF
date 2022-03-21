/*! \file utilities.cxx
 *  \brief this file contains an assortment of utilities
 */

#include "ioutils.h"
#include "logging.h"
#include "stf.h"

namespace vr {

std::string basename(const std::string &filename)
{
    std::string basenamed(filename);
    auto last_sep = basenamed.rfind('/');
    if (last_sep != std::string::npos) {
        basenamed = basenamed.substr(last_sep + 1);
    }
    return basenamed;
}

} // namespace vr

int CompareInt(const void *p1, const void *p2) {
      Int_t val1 = *(Int_t*)p1;
      Int_t val2 = *(Int_t*)p2;
      return (val1 - val2);
}

/// get the memory use looking at the task
namespace vr {
    memory_usage get_memory_usage() {
        memory_usage usage;

        const char *stat_file = "/proc/self/status";
        std::ifstream f(stat_file);
        if (!f.is_open()) {
            LOG(warning) << "Couldn't open " << stat_file << " for memory usage reading";
        }
        for (std::string line; std::getline(f, line); ) {
            auto start = line.find("VmSize:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 7));
                is >> usage.vm.current;
                continue;
            }
            start = line.find("VmPeak:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 7));
                is >> usage.vm.peak;
                continue;
            }
            start = line.find("VmRSS:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 6));
                is >> usage.rss.current;
                continue;
            }
            start = line.find("VmHWM:");
            if (start != std::string::npos) {
                std::istringstream is(line.substr(start + 6));
                is >> usage.rss.peak;
                continue;
            }
        }

        // all values above are in kB
        usage.vm.current *= 1024;
        usage.vm.peak *= 1024;
        usage.rss.current *= 1024;
        usage.rss.peak *= 1024;
        return usage;
    }
} // namespace vr

std::string GetMemUsage(const std::string &function)
{
    auto memory_usage = vr::get_memory_usage();
    std::ostringstream memory_report;
    auto append_memory_stats = [&memory_report](const char *name, const vr::memory_stats &stats) {
        memory_report << name << " current/peak: " << vr::memory_amount(stats.current) << " / " << vr::memory_amount(stats.peak);
    };
    memory_report << "Memory report @ " << function << ": ";
    append_memory_stats("VM", memory_usage.vm);
    memory_report << "; ";
    append_memory_stats("RSS", memory_usage.rss);
    return memory_report.str();
}

/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
#ifdef __APPLE__

static inline void
CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

static inline void
CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

static inline int
CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set)
{
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (int i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}
#endif

void cpuset_to_cstr(cpu_set_t *mask, char *str)
{
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d ", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d ", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d ", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  ptr = nullptr;
}

void report_binding()
{
    // if there is no MPI do not report any binding
#if !defined(USEMPI) && !defined(USEOPENMP)
    return; 
#endif
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#else
    int ThisTask=0, NProcs=1;
#endif
    string binding_report = "Core binding \n ";
    cpu_set_t coremask;
    char clbuf[7 * CPU_SETSIZE], hnbuf[64];
    memset(clbuf, 0, sizeof(clbuf));
    memset(hnbuf, 0, sizeof(hnbuf));
    (void)gethostname(hnbuf, sizeof(hnbuf));
#ifdef USEOPENMP
    #pragma omp parallel shared (binding_report) private(coremask, clbuf) 
#endif
    {
        string result;
        (void)sched_getaffinity(0, sizeof(coremask), &coremask);
        cpuset_to_cstr(&coremask, clbuf);
        result = "\t On node " + string(hnbuf) + " : ";
#ifdef USEMPI 
        result += "MPI Rank " + to_string(ThisTask) + " : ";
#endif
#ifdef USEOPENMP
        auto thread = omp_get_thread_num();
        result +=" OMP Thread " + to_string(thread) + " : ";
#endif
        result += " Core affinity = " + string(clbuf) + " \n ";
#ifdef USEOPENMP 
        #pragma omp critical 
#endif 
        {
            binding_report +=result;

        }
    }
    LOG(info)<<binding_report;
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}


#ifdef NOMASS
void VR_NOMASS(){};
#endif
#ifdef GASON
void VR_GASON(){};
#endif
#ifdef STARON
void VR_STARON(){};
#endif
#ifdef BHON
void VR_BHON(){};
#endif
#ifdef USEMPI
void VR_MPION(){};
#endif
#ifdef USEOPENMP
void VR_OPENMPON(){};
#endif
#ifdef HIGHRES
void VR_ZOOMSIMON(){};
#endif
#ifdef USEHDF
void VR_HDFON(){};
#ifdef USEPARALLELHDF
void VR_PARALLELHDFON(){};
#endif
#endif
