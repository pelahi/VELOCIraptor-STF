/*! \file ompvar.h
 *  \brief this file declares variables used for openmp compilation.
 */

#ifndef OMPVAR_H
#define OMPVAR_H

#ifdef USEOPENMP
#include <omp.h>
#endif

///\name Include for NBodyFramework library.
//@{
///nbody code
#include <NBody.h>
///Math code
#include <NBodyMath.h>
///Binary KD-Tree code
#include <KDTree.h>
//@}

using namespace std;
using namespace Math;
using namespace NBody;

/// \defgroup OMPLIMS For determining whether loop contains enough for openm to be worthwhile.
//@{
#define ompsplitsubsearchnum 10000000
#define ompsubsearchnum 10000
#define ompsearchnum 50000
#define ompunbindnum 1000
#define ompperiodnum 1000000
#define omppropnum 50000
#define ompfofsearchnum 2000000
#define ompsortsize 1000000
//@}

#ifdef USEOPENMP 

///structure to store relevant info for searching openmp domains
struct OMP_Domain {
    Int_t ncount, noffset, numgroups;
    Double_t bnd[3][2];
    vector<int> neighbour;
};

///structure relevant for exchanging particle information between openmp domains
struct OMP_ImportInfo {
    Int_t index, pfof;
    int task;
};

/// Useful class for managing threads and gpus
class VROMPThreadPool
{
    public :
        /// allow for parallel omp
        bool allowomp, allowoffload;
        /// total number of levels 
        unsigned int nlevels;
        /// current omp levle 
        unsigned int ncurlevel; 
        ///total number of threads 
        unsigned int nthreads;
        ///number of active idle threads that can be given tasks to do 
        unsigned int nactivethreads;
        ///total number of gpus 
        unsigned int ngpus;
        ///number of active idle gpus that can be given tasks to do 
        unsigned int nactivegpus;
        ///vector storing ids of active idle threads and gpus, 
        vector<unsigned int> activethreadids, activegpuids;
        ///vector storing ids of active threads and gpus, 
        vector<unsigned int> idlethreadids, idlegpuids;
        ///hiearchy level of pool. Pools of the same hierarchy cannot be merged 
        ///then a pool is split, a child power of a higher hierarchy level is produced 
        unsigned int hierarchylevel;

        /// \name Constructors 
        //@{
        VROMPThreadPool();
        VROMPThreadPool(const VROMPThreadPool &) = default;
        VROMPThreadPool(VROMPThreadPool &&) = default;
        /// build pool based on input VROMPThreadPool pointer or 
        /// system defaults if pointer to thread pool null 
        VROMPThreadPool(VROMPThreadPool *);
        VROMPThreadPool& operator=(const VROMPThreadPool&) = default;
        VROMPThreadPool& operator=(VROMPThreadPool&&) = default;
        //@}

        /// \name Thread Management Routines
        //@{
        /// Init the thread and gpu pool 
        void Init();
        /// split the resources and produce new pool
        VROMPThreadPool Split();
        VROMPThreadPool Split(unsigned int nthreadsplit, unsigned int ngpusplit);
        /// merge resources from a 
        void Merge(VROMPThreadPool &);
        /// move thread from idle pool to active pool 
        void ActivateThread(unsigned int nactivate=1);
        ///  deactivate a thread and move it to idle pool  
        void DeactivateThread(unsigned int ndeactivate=1);
        /// move GPU from idle pool to active pool 
        void ActivateGPU(unsigned int nactivate=1);
        ///  deactivate a GPU and move it to idle pool  
        void DeactivateGPU(unsigned int ndeactivate=1);

        /// close pool 
        void Close();


        //@}

        void Print();

};

#endif
#endif
