/*! \file mpivar.h
 *  \brief this file declares variables used for mpi compilation.
 */

#ifndef MPIVAR_H
#define MPIVAR_H

#include <mpi.h>

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

/// \name for mpi tasks and domain construction
//@{
extern int ThisTask, NProcs;
//@}

#endif
