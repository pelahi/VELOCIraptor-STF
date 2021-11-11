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

void InitMemUsageLog(Options &opt){
#ifndef USEMPI
    int ThisTask=0;
#endif
    if (!opt.memuse_log) return;
    ofstream Fmem;
    char buffer[2500];
    sprintf(buffer,"%s.memlog.%d",opt.outname,ThisTask);
    Fmem.open(buffer);
    Fmem<<"Memory Log"<<endl;
    Fmem.close();
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
