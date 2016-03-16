/*! \file allvars.cxx
 *  \brief provides instances of global variables
 */
#include "allvars.h"

#ifdef USEMPI
/// \name For MPI 
//@{
int ThisTask,NProcs;
int NSnap,StartSnap,EndSnap;
//@}

#endif
