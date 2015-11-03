#ifndef DENSITY_H
#define DENSITY_H

#include <NBody.h>
#include <NBodyMath.h>
namespace NBody
{
    /// Bin system in shperical shells, and compute density in each shell.
    /// Bin spacing is normal.    
    void DensityBinNorm(const System &S, int NumBins, Double_t *RBin, 
                        Double_t *DBin);
    
    /// Bin system in shperical shells, and compute density in each shell.
    /// Bin spacing is logarithmic.
    void DensityBinLog(const System &S, int NumBins, Double_t *RBin, 
                       Double_t *DBin);
    
    /// Bin system in spherical shells that have the same number of particles
    /// in each.
    void DensityBinEqualMass(System &S, int NumBins, Double_t *RBin, 
                             Double_t *DBin);
    
    /// Calculate the density of the system in elliptical bins.
    /// q and s are axis ratios of the system.
    void DensityBinLogEllip(const System &S, int NumBins, Double_t *RBin, 
                            Double_t *DBin, Double_t q, Double_t s);
    
    /// Calculate the density of each particle by smoothing over
    /// nearest neighbours.  Uses a KD-tree to locate nearest neighbours,
    /// and a smoothing kernel to calculate density.  Density is stored
    /// per particle (field rho).
    void DensitySmooth(System &S, int BucketSize = 16, int NumNearest = 64);

    /// Calculate the density of each particle by smoothing all neighbours within a radius R
    void DensitySmoothBall(System &S, Double_t R);
}  

#endif
