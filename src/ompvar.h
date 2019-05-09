/*! \file ompvar.h
 *  \brief this file declares variables used for openmp compilation.
 */

#ifndef OMPVAR_H
#define OMPVAR_H

#include <omp.h>

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
#define ompsearchnum 50000
#define ompunbindnum 1000
#define ompperiodnum 100000
#define omppropnum 50000
#define ompfofsearchnum 2000000
//@}

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
