/*! \file System.h
 *  \brief header file for the \ref NBody::System Class

*/

#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <NBodyMath.h>
#include <Particle.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Math;
namespace NBody
{

/*! 
    \class NBody::System
    \brief A system of particles
    
    This class represents a system of particles. By default it builds an array of \ref NBody::Particle, but also can allocate 
    memory for \ref NBody::GasParticle and \ref NBody::StarParticle. These classes are defined in \ref Particle.h. 
*/

  class System
  {
    protected:
    /// \name Current number of particles in system (also gas, star)
    //@{
    Int_t numparts;
    Int_t ndark, ngas,nstar;
    //@}
    /// Time stamp of system
    Double_t time;
    /// \name An array of particles (mass, x, y, z, vx, vy, vz,id,type, rho), gas particles, and star particles
    //@{
    Particle *particle;
    GasParticle *gparticle;
    StarParticle *sparticle;
    //@}
    ///Period of box (if zero not periodic)
    Coordinate period;

    public:

    /// \name Constructors & Destructors
    //@{
    System(Int_t num = 1, Double_t t = 0.0, Double_t *period=NULL, Int_t ng=0, Int_t ns=0);
    System(Int_t num, Particle *p, Double_t t = 0.0, Double_t *period=NULL, Int_t ng=0, Int_t ns=0, GasParticle *gp=NULL, StarParticle *sp=NULL);
    System(const System &s);
    ~System() { delete [] particle; if (ngas>0) delete[] gparticle; if (nstar) delete[] sparticle;}
    //@}

    /// \name Get & Set methods
    //@{
    Int_t GetNumParts() const { return numparts; }
    Int_t GetNumDark() const { return ndark; }
    Int_t GetNumGas() const { return ngas; }
    Int_t GetNumStar() const { return nstar; }

    Double_t GetTime() const { return time; }
    void SetTime(Double_t t) { time=t; }

    Particle *Parts() const { return particle; }
    Particle &Part(Int_t i) const { return particle[i]; }

    Coordinate GetPeriod() const {return period;}
    void SetPeriod(Coordinate p){period=p;}
    void SetPeriod(Coord p){for (int i=0;i<3;i++)period[i]=p.pos[i];}
    void SetPeriod(Double_t *p){for (int i=0;i<3;i++)period[i]=p[i];}
    //@}

    //    OTHER FUNCTIONS

    /// \name Properties of the system
    ///in these functions bool variable indicates whether system has already been adjust to CM 
    ///frame needed for some of the functions.
    //@{

    Double_t TotalMass() const;
    ///returns maximum x,y, or z coord, ie length of box needed to
    ///contain most distant particle
    Double_t MaxLength() const;
    ///returns maximum radius for sphere needed to
    ///contain most radially distant particle from CM
    Double_t MaxRadius(bool cmframe=true) const;
    ///as above but returns minimum
    Double_t MinRadius(bool cmframe=true) const;
    ///as above but returns limits (min, average, max)
    Double_t *RadialLimits(bool cmframe=true) const;
    ///returns the average distance from the CM
    Double_t AverageRadius(bool cmframe=true) const;
    ///returns the average speed in the CM frame
    Double_t AverageSpeed(bool cmframe=true) const;
    ///returns the average velocity in the CM frame
    Coordinate AverageVelocity(bool cmframe=true) const;
    ///returns the average radial velocity in the CM frame
    Double_t AverageRadialVelocity(bool cmframe=true) const;
    ///returns the maximum radial velocity in the CM frame
    Double_t MaxRadialVelocity(bool cmframe=true) const;
    ///returns the maximum circular velocity in the CM frame
    Double_t MaxCircularVelocity(bool cmframe=true) const;
    ///returns the average density of the system
    Double_t AverageDensity() const;
    ///returns the density in some radius R from some point O 
    Double_t RegionDensity(Double_t R, Double_t x, Double_t y, Double_t z) const;
    ///Angular Momentum (one with argument is angular momentum about a point)
    ///one without is angular momentum about center of mass
    Coord AngularMomentum() const;
    Coord AngularMomentum(Coord) const;
    Coord AngularMomentum(Coordinate) const;

    //@}

    /// \name Sorting Functions
    //@{
    
    ///Sort particles in system by radius
    void SortByID();
    ///Sort particles in system by radius
    void SortByRadius();
    //void SortByRadius();
    //sorts from index start  n particles.
    //void SortByRadius(int start, int n);
    ///Sort particles in system by distance to some point
    void SortByDistance(Coordinate c);
    ///sorts from index start  n particles.
    void SortByDistance(Coordinate c, int start, int n);
    ///Sort system by energy (most bound to least bound)
    void SortByEnergy(Double_t G, Double_t eps);
    ///sorts from index start  n particles.
    void SortByEnergy(Double_t G, Double_t eps, int start, int n);
    ///Sort particles in system by Density
    void SortByDensity();
    ///sorts from index start  n particles.
    void SortByDensity(int start, int n);
    
    //@}
    
    /// \name CM and reference frame functions
    //@{
    
    ///returns center of mass
    Coord CM(Double_t tolerance=0.0) const;
    ///Subtract centre of mass coords
    void AdjustForCM(Double_t tolerance=0.0);
    ///returns centre of velocity
    Coord CMVel(Double_t tolerance=0.0) const;
    ///Subtract centre of mass velocity
    void AdjustForCMVel(Double_t tolerance=0.0);
    /// Subtract coordinate, essentially specify new origin (non-relatvistic)
    void AdjustPosition(Coord);
    /// Subtract coordinate, essentially specify new origin (non-relatvistic)
    void AdjustPosition(Coordinate);
    /// Subtract velocity, essentially specify new frame (non-relativistic)
    void AdjustVelocity(Coord);
    /// Subtract velocity, essentially specify new frame (non-relativistic)
    void AdjustVelocity(Coordinate);

    //@}

    /// \name Energy functions
    ///Note that this KE has units of M*L^2/T^2
    //@{

    ///Calculates the kinetic energy of  particle i about the CM
    Double_t KineticEnergy(Int_t i) const;
    ///Calculates the kinetic energy of the particles about the CM
    Double_t KineticEnergy() const;
    ///Calculates the kinetic energy of  particle i in frame with a velocity
    Double_t KineticEnergy(Coord, Int_t ) const;
    ///Calculates the kinetic energy of the particles in frame with a velocity
    Double_t KineticEnergy(Coord) const;

    ///Calculates the potential energy of particle between the ith and jth particle of the System
    Double_t PotentialEnergy(Int_t i, Int_t j, Double_t eps=0.0) const;
    ///Calculates the potential energy of particle i in the well of the System
    Double_t PotentialEnergy(Int_t i, Double_t eps=0.0) const;
    ///Calculates the potential energy of the particles in the gravtiational
    ///well of the System
    Double_t PotentialEnergy(Double_t eps=0.0, bool sphericalapprox=false) const;

    //@}

    /// \name Manipulate particles in system
    //@{
    
    ///Finds a particle and returns the indice.
    Int_t FindParticle(const Particle &p) const;
    
    /// Add/remove particle to/from system, if possible.
    void AddParticle(Particle part);
    void RemoveParticle(Particle part);
    void RemoveParticle(Int_t i);
    ///Adds a system based on a vector between origins.
    ///o must be position of origin of system to be added in coords of this.
    void AddSystem(const System &s, const Double_t* o);
    ///vo must be the velocity of the origin of the system in coords of this.
    void AddSystem(const System &s, const Double_t* o, const Double_t* vo);
    /// Remove all particles beyond a certain radius
    void ExtractSphere(Double_t radius);
    //do not plan at the moment to implements this. Very tricky and also
    //not likely to be needed.
    //void RemoveSystem(const System &s);
    
    //@}


    /// \name Get & Set methods
    //@{
    System &operator=(const System &);
    bool operator==(const System &) const;

    //To simplisitic should not allow. Force user to use AddSystem;
    //System &operator+(const System &);
    //System &operator-(const System &);
    //Force users for now to use add/remove particle functions.
    //Also think code below this is incorrect anyway.
    //System &operator+(const Particle &p) { AddParticle(p); }
    //System &operator-(const Particle&p) { RemoveParticle(p); }

    Particle &operator[](Int_t i);
    Particle operator[](Int_t i) const;

    /// Inserter: prints system as ascii data
    friend ostream &operator << (ostream &outs, const System &s);
    //@}

    // Extractor:  get system from ascii data
    // Warning: very slow
    //friend istream &operator >> (istream &ins, System &s);

    /// \name IO functions
    //@{
    
    //Read system from binary data
    // ... to a file stream ...
    //void Read(ifstream &F);
    //Output system as binary data
    // ... to a file stream ...
    void Write(ofstream &F)const;
    // ... to a file ...
    void Write(char *Filename)const;
    // ... or to the screen.  This is c-style io.
    void Write(FILE *stream = stdout)const;
    
    //@}

  };

  /*
      End of System Class
  */
    /*
        Sub Classes: Gas and Star systems ///not really implemented yet
    class GasSystem : public System
    {
        private:
            Int_t ngas;
            GasParticle *gparticle;
        public:

        Double_t TotalGasMass() const;
        Double_t AverageTemp() const;
        Double_t MaxTemp() const;
    };
    */
}

#endif // SYSTEM_H
