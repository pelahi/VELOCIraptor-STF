#ifndef ENERGY_H
#define ENERGY_H

#include <NBody.h>

namespace NBody
{
    /// The type of orbits in the system.
    enum Orbit { Radial, Isotropic };
        
    /// Calculate the gravitational potential energy by direct summation.
    /// Don't use this for more than ~10000 particles!
    void CalcPotDirect(System &S, Double_t eps = 0.0);
    
    /// Calculate the gravitational potential energy by assuming the
    /// system is spherically symmetric.  The system MUST be sorted
    /// by radius before this function is called.
    void CalcPotShell(System &S, Double_t eps = 0.0);
    
    /// Calculate the energy of each particle, E = v^2/2 + Phi.
    /// Phi should be calculate and stored in the particle data.
    void CalcEnergy(const System &S, Double_t *E, Double_t &Etot);
    
    /// Calculate the differential energy distribution
    void DiffEnergyDistNorm(const System &S, const Double_t *E, Double_t *Ebin,
                            Double_t *Mbin, int NumBins);
    
    /// Calculate the differential energy distribution using bins of equal
    /// mass.
    void DiffEnergyDistEqual(System &S, const Double_t *E, Double_t *Ebin,
                             Double_t *Mbin, int NumBins);
    
    /// Calculate the number density N(E, j).  Nbin is an array
    /// of length NumBins * NumBins.
    void NumberDensity(const System &S, const Double_t *E, Double_t *Ebin, 
                       Double_t *Jbin, Double_t *Nbin, int NumBins);
    
    /// Calculate the density of states, g(E)
    void DenStates(const System &S, Double_t *g, const Double_t *Ebin, 
                   int NumBins, Orbit orbit);

}  

#endif
