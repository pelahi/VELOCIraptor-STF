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
#endif

/*
 * OMP_Threadpool_datastructure
 * - coupled to OMP_Threadpool
 * - structure to store threadpool information
 * - captures the following:
 *      thread state active(True)/inactive(false)
 *      total number of particles
 */
struct OMP_Threadpool_datastructure {
    bool active;
    Int_t particle_total;
};

/*
 * OMP_Threadpool
 * - contains routines used for optimised parallel computing/processing
 *      throttle_up - increase number of active threads
 *      throttle_down - decrease number of active threads
 * - captures the following:
 *      specific threadpool information stored to OMP_Threadpool_datastructure
 *      total thread size
 *      number of active threads
 */
struct OMP_Threadpool {
    // Total thread available
    // Active threads
    // Unsigned long long for number of particles
    // omp thread class - 
    // omp retrieve active threads
    // specify offset of thread

    /* Private block */
    private:
    map<int, OMP_Threadpool_datastructure> details;
    int thread_count, active_threads;

    /* Public block */
    public:
    void throttle_up(){
        // TODO: add the following logic:
        //      OMP_Threadpool_datastructure update state of thread from inactive to active
        active_threads += 1;
    }

    void throttle_down(){
        // TODO: add the following logic:
        //      OMP_Threadpool_datastructure update state of thread from active to inactive
        active_threads -= 1;
    }

    void set_total_threads(int thread_count){
        thread_count = thread_count;
    }

    int get_total_threads(){
        return thread_count;
    }

    int get_active_threads(){
        return active_threads;
    }

    /* Confirmation Functions */
    void print_total_threads(){
        printf("Total Threads: %d\n", thread_count);
    }

    void print_active_threads(){
        printf("Active Threads: %d\n", active_threads);
    }
};
