/*! \file omproutines.cxx
 *  \brief this file contains routines used with OpenMP compilation.

    OpenMP routines pertaining to complex OpenMP calls.
 */

#ifdef USEOPENMP

//-- For MPI

#include "stf.h"

/// \name routines which check to see if some search region overlaps with local mpi domain
//@{
///search if some region is in the local mpi domain
int OpenMPSearchForOverlap(Double_t xsearch[3][2], Double_t bnd[3][2], Double_t period){
    Double_t xsearchp[3][2];
    if (!((bnd[0][1] < xsearch[0][0]) || (bnd[0][0] > xsearch[0][1]) ||
        (bnd[1][1] < xsearch[1][0]) || (bnd[1][0] > xsearch[1][1]) ||
        (bnd[2][1] < xsearch[2][0]) || (bnd[2][0] > xsearch[2][1])))
            return 1;
    else {
        if (period==0) return 0;
        else {
            for (int j=0;j<3;j++) {xsearchp[j][0]=xsearch[j][0];xsearchp[j][1]=xsearch[j][1];}
            for (int j=0;j<3;j++) {
                if (!((bnd[j][1] < xsearch[j][0]+period) || (bnd[j][0] > xsearch[j][1]+period))) {xsearchp[j][0]+=period;xsearchp[j][1]+=period;}
                else if (!((bnd[j][1] < xsearch[j][0]-period) || (bnd[j][0] > xsearch[j][1]-period))) {xsearchp[j][0]-=period;xsearchp[j][1]-=period;}
            }
            if (!((bnd[0][1] < xsearchp[0][0]) || (bnd[0][0] > xsearchp[0][1]) ||
            (bnd[1][1] < xsearchp[1][0]) || (bnd[1][0] > xsearchp[1][1]) ||
            (bnd[2][1] < xsearchp[2][0]) || (bnd[2][0] > xsearchp[2][1])))
                return 1;
            else return 0;
        }
    }
}

/// Determine if a particle needs to be exported to another mpi domain based on a physical search radius
int OpenMPSearchForOverlap(Particle &Part, Double_t bnd[3][2], Double_t rdist, Double_t period){
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=Part.GetPosition(k)-rdist;xsearch[k][1]=Part.GetPosition(k)+rdist;}
    return OpenMPSearchForOverlap(xsearch,bnd,period);
}

//@}

#endif
