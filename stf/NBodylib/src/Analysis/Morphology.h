#ifndef MORPHOLOGY_H
#define MORPHOLOGY_H

#include <NBodyMath.h>
#include <NBody.h>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Math;
namespace NBody
{
    /// Calculate the global morphology of the system.  See, e.g., Katz
    /// (1991) or Dubinski and Carlberg (1991).  q and s are the axis ratios 
    /// of the system.  Error is relative error for convergence.
	/// assumes equal masses
    void GetGlobalMorphology(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, int verbose=0);
    //also returns the eigenvectors.
    void GetGlobalMorphology(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, Matrix& Eigenvectors, int verbose=0);
	/// does not assume equal mass (longer computation time and masses for convergence should be order unity)
    void GetGlobalMorphologyWithMass(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, int mdenom=0, int verbose=0);
    //also returns the eigenvectors.
    void GetGlobalMorphologyWithMass(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, Matrix& Eigenvectors, int mdenom=0, int verbose=0);

	/// Same as above with different interface
    void GetGlobalMorphology(System& S, Double_t& q, Double_t& s, Double_t Error, int verbose=0);
    void GetGlobalMorphology(System& S, Double_t& q, Double_t& s, Double_t Error, Matrix& Eigenvectors, int verbose=0);
    void GetGlobalMorphologyWithMass(System& S, Double_t& q, Double_t& s, Double_t Error, int mdenom=0, int verbose=0);
    void GetGlobalMorphologyWithMass(System& S, Double_t& q, Double_t& s, Double_t Error, Matrix& Eigenvectors, int mdenom=0, int verbose=0);

	/// general calculation of interia tensor
    void GetInertiaTensor(const Int_t nbodies, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigvec);
    void GetInertiaTensor(const Int_t nbodies, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigvec, Matrix &I);
	/// Same as above with different interface
    void GetInertiaTensor(System& S, Double_t &a, Double_t &b, Double_t &c, Matrix& eigvec);
    void GetInertiaTensor(System& S, Double_t &a, Double_t &b, Double_t &c, Matrix& eigvec, Matrix &I);
//    void GetInertiaTensor(System& S, Double_t &q, Double_t &s, Matrix& eigvec);
}

#endif // MORPHOLOGY_H
