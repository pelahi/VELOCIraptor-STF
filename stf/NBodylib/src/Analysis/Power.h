#ifndef POWER_H
#define POWER_H

#include <NBody.h>

namespace NBody
{
    /// Calculate the density on a cubic grid.  This works best
    /// if the system is cubic, of course.  The 1-D array rho has
    /// size N^3, and the ijk'th element is accessed as:
    /// rho_ijk = rho[k + N*(j + N*i)].  Uses nearest-grid-point.
    /// Also computes the box length for you.  Warning:  this currently
    /// expects a box that goes from 0 ... L.
    void DensityGridNGP(System &S, Double_t *rho, int Ng, Double_t &L);
    
    /// Calculate the power spectrum of a system.  rho is an array
    /// of length [Ng * Ng * Ng] containing gridded density information.
    void CalcPowSpectrum(Double_t *rho, Double_t *p, int Ng, Double_t BoxLength);
}  

#endif
