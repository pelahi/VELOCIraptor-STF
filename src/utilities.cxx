/*! \file utilities.cxx
 *  \brief this file contains an assortment of utilities
 */

#include "stf.h"

/// Time functions
//@{
int GetMilliCount()
{
  // Something like GetTickCount but portable
  // It rolls over every ~ 12.1 days (0x100000/24/60/60)
  // Use GetMilliSpan to correct for rollover
  timeb tb;
  ftime( &tb );
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}

int GetMilliSpan( int nTimeStart )
{
  int nSpan = GetMilliCount() - nTimeStart;
  if ( nSpan < 0 )
    nSpan += 0x100000 * 1000;
  return nSpan;
}

//@}

int CompareInt(const void *p1, const void *p2) {
      Int_t val1 = *(Int_t*)p1;
      Int_t val2 = *(Int_t*)p2;
      return (val1 - val2);
}

/// get the memory use looking at the task
void GetMemUsage(Options &opt, string funcname, bool printreport){
#ifndef USEMPI
    int ThisTask=0;
#endif
    // if (!opt.memuse_log) return;
    ifstream f;
    size_t pos1, pos2;
    unsigned long long size, resident, shared, text, data, library, dirty, peak;
    map<string,float> memuse;
    bool iflag = true;
    char buffer[2500];
    string memreport, temp, delimiter1("VmPeak:"), delimiter2("kB");
    ofstream Fmem;
    // Open the file storing memory information associated with the process;
    f.open("/proc/self/statm");
    if (f.is_open()) {
        f >> size >> resident >> shared >> text >> library >> data >>dirty;
        // nscan = fscanf(file, "%lld %lld %lld %lld %lld %lld %lld",
        //     &size, &resident, &shared, &text, &library, &data, &dirty);
        f.close();
    }
    else {
        iflag = false;
    }
    f.open("/proc/self/status");
    if (f.is_open()) {
        while (f.getline(buffer, 2500)){
            temp = string(buffer);
            if ((pos1 = temp.find(delimiter1)) != string::npos) {
                pos2 = temp.find(delimiter2);
                temp = temp.substr(pos1+delimiter1.size(), pos2);
                peak = stol(temp)*1024;
                break;
            }
        }
        f.close();
    }
    else {
        iflag = false;
    }

    memreport += string("Memory report, func = ")+funcname+string(" task = ")+to_string(ThisTask)+ string(" : ");

    //having scanned data for memory footprint in pages, report
    //memory footprint in GB
    if (iflag) {
        // Convert pages into bytes. Usually 4096, but could be 512 on some
        // systems so take care in conversion to KB. */
        uint64_t sz = sysconf(_SC_PAGESIZE);
        size *= sz ;
        resident *= sz ;
        shared *= sz ;
        text *= sz ;
        library *= sz ;
        data *=  sz ;
        dirty *= sz ;
        // float bytestoGB = 1.0/(1024.0*1024.*1024.);
        float bytestoGB = 1.0;
        if (opt.memuse_peak < peak) opt.memuse_peak = peak;
        opt.memuse_nsamples++;
        opt.memuse_ave += size;
        memuse["Size"] = size*bytestoGB;
        memuse["Resident"] = resident*bytestoGB;
        memuse["Shared"] = shared*bytestoGB;
        memuse["Text"] = text*bytestoGB;
        memuse["Library"] = library*bytestoGB;
        memuse["Data"] = data*bytestoGB;
        memuse["Dirty"] = dirty*bytestoGB;
        memuse["Peak"] = opt.memuse_peak*bytestoGB;
        memuse["Average"] = opt.memuse_ave/(float)opt.memuse_nsamples*bytestoGB;

        for (map<string,float>::iterator it=memuse.begin(); it!=memuse.end(); ++it) {
            memreport += it->first + string(" = ") + to_string(it->second) + string(" GB, ");
        }
    }
    else{
        memreport+= string(" unable to open or scane system file storing memory use");
    }
    // sprintf(buffer,"%s.memlog.%d",opt.outname,ThisTask);
    // Fmem.open(buffer,ios::app);
    // Fmem<<memreport<<endl;
    // Fmem.close();
    if (printreport) cout<<memreport<<endl;
}

/// get the memory use looking at the task
void GetMemUsage(string funcname, bool printreport){
#ifndef USEMPI
    int ThisTask=0;
#endif
    ifstream f;
    size_t pos1, pos2;
    unsigned long long size, resident, shared, text, data, library, dirty, peak;
    map<string,float> memuse;
    bool iflag = true;
    char buffer[2500];
    string memreport, temp, delimiter1("VmPeak:"), delimiter2("kB");
    ofstream Fmem;
    // Open the file storing memory information associated with the process;
    f.open("/proc/self/statm");
    if (f.is_open()) {
        f >> size >> resident >> shared >> text >> library >> data >>dirty;
        // nscan = fscanf(file, "%lld %lld %lld %lld %lld %lld %lld",
        //     &size, &resident, &shared, &text, &library, &data, &dirty);
        f.close();
    }
    else {
        iflag = false;
    }
    f.open("/proc/self/status");
    if (f.is_open()) {
        while (f.getline(buffer, 2500)){
            temp = string(buffer);
            if ((pos1 = temp.find(delimiter1)) != string::npos) {
                pos2 = temp.find(delimiter2);
                temp = temp.substr(pos1+delimiter1.size(), pos2);
                peak = stol(temp)*1024;
                break;
            }
        }
        f.close();
    }
    else {
        iflag = false;
    }

    memreport += string("Memory report, func = ")+funcname+string(" task = ")+to_string(ThisTask)+ string(" : ");

    //having scanned data for memory footprint in pages, report
    //memory footprint in GB
    if (iflag) {
        // Convert pages into bytes. Usually 4096, but could be 512 on some
        // systems so take care in conversion to KB. */
        uint64_t sz = sysconf(_SC_PAGESIZE);
        size *= sz ;
        resident *= sz ;
        shared *= sz ;
        text *= sz ;
        library *= sz ;
        data *=  sz ;
        dirty *= sz ;
        // float bytestoGB = 1.0/(1024.0*1024.*1024.);
        float bytestoGB = 1.0;
        // if (opt.memuse_peak < peak) opt.memuse_peak = peak;
        // opt.memuse_nsamples++;
        // opt.memuse_ave += size;
        memuse["Size"] = size*bytestoGB;
        memuse["Resident"] = resident*bytestoGB;
        memuse["Shared"] = shared*bytestoGB;
        // memuse["Text"] = text*bytestoGB;
        // memuse["Library"] = library*bytestoGB;
        memuse["Data"] = data*bytestoGB;
        // memuse["Dirty"] = dirty*bytestoGB;
        memuse["Peak"] = peak*bytestoGB;
        // memuse["Average"] = opt.memuse_ave/(float)opt.memuse_nsamples*bytestoGB;

        for (map<string,float>::iterator it=memuse.begin(); it!=memuse.end(); ++it) {
            memreport += it->first + string(" = ") + to_string(it->second) + string(" GB, ");
        }
    }
    else{
        memreport+= string(" unable to open or scane system file storing memory use");
    }
    if (printreport) cout<<memreport<<endl;
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

double MyGetTime(){
#ifdef USEOPENMP
    return omp_get_wtime();
#else
    return (clock() /( (double)CLOCKS_PER_SEC));
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
