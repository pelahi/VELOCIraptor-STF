#ifndef INITCOND_H
#define INITCOND_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

//nbody code
#include <NBody.h>
#include <NBodyMath.h>

using namespace Math;
using namespace NBody;
using namespace std;

namespace NBody
{
    //simple initial conditions.

    /// Sets up the system as a cubic grid of particles, with
    /// zero velocities.  N here is the number of particles per side.
    /// origin is at corner.
    System* CubicGrid(int N, Double_t BoxLength, Double_t Mtotal,
                      Double_t time = 0.0);
    
    /// Sets up the system as a uniform sphere, with
    /// zero velocities.  N here is the initial number of particles.
    /// sphere is centered.
    System* UniformSphere(int N, Double_t BoxLength, Double_t Mtotal,
                      Double_t time = 0.0);
    /// Calculate Gaussian displacement vectors from the power
    /// spectrum Pow.  Vectors Si are arrays such that
    /// Si(x,y,z) = Si[z + N * (y + N * x)]
    //void CalcGaussDisp(int N, Double_t BoxLength, Double_t *Sx, Double_t *Sy,
    //                   Double_t *Sz, Cosmology &cosmo, unsigned int Seed = 0);
    
    /// Displace particles, initially on a uniform grid, using the
    /// Zeldovich approximation.
    //void ApplyZelApprox(System &s, const Double_t *Sx, const Double_t *Sy,
    //                    const Double_t *Sz, Double_t BoxLength, Cosmology &comso,
    //                    bool GadgetVpec = false);
    
    /// Create a Plummer system.  N is the total number of particles.
    /// Useful for testing things.  To come.
    //System* Plummer(int N, Double_t t = 0.0);
}  

#endif
