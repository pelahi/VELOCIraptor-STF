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
#endif
