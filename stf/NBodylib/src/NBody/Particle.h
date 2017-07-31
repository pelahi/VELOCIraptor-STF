/*! \file Particle.h
 *  \brief header file for the \ref NBody::Particle and derived classes

*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <NBodyMath.h>


#ifdef USEBOOSTMPI
#include <boost/mpi.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#endif

using namespace std;
using namespace Math;
namespace NBody
{


/// \name Particle Memory Flags
/// \brief These are define flags that alter the memory allocation of particles
//@{
///use long ints to store particle ids
///if NOMASS is defined, don't actually store a mass
//#define NOMASS
///this is the mass value that is returned if no value is actually stored
#define MASSVAL 1.0
///use floats to store positions and velocities regardless of the precision used elsewhere
#ifdef LOWPRECISIONPOS
typedef float DoublePos_t;
#else
typedef Double_t DoublePos_t;
#endif
//@}

/*!
    \class NBody::Particle
    \brief A simple n-body particle class.

    The class is for a simple n-body particle with mass, position, velocity, id, type and density information.
    Note that here particles can be compiled to store mass or assume all particles have the same mass given by a
    MASSVAL definition in \ref Particle.h by defining the flag NOMASS. Particles can also be compiled to have long ids
    stored using an unsigned int or unsigned long it.
*/
    class Particle
    {
        protected:
        ///mass
#ifndef NOMASS
        Double_t mass;
#endif
#ifdef LOWPRECISIONPOS
        /// x, y, z
        DoublePos_t position[3];
        /// vx, vy, vz
        DoublePos_t velocity[3];
#else
        /// x, y, z
        Double_t position[3];
        /// vx, vy, vz
        Double_t velocity[3];
#endif
        ///\name identifiers
        //@{
        //extra particle identifier, useful as it is not changed by KDTree
#ifdef LONGPIDS
        UInt_t pid;
#else
        Int_t pid;
#endif
#ifdef LONGIDS
        UInt_t id;
#else
        Int_t id;
#endif
        //int gid;
        int type;
        //@}
        ///density
        Double_t rho;
        ///potential
        Double_t phi;
        ///For hydrodynamical quantities
        //@{
        ///if gas flag is set, then particle also can have sph based quantities. if star flag is set,
        ///then particle also can have star properties such as age and metallicity
#ifdef GASON
        ///self energy
        DoublePos_t u;
        ///sph based density
        DoublePos_t sphden;
#endif
#ifdef STARON
        ///stellar age
        DoublePos_t tage;
#endif
#if defined (GASON) && (STARON)
        ///metallicity
        DoublePos_t zmet;
        ///star formation rate of gas
        DoublePos_t sfr;
#endif
#if defined(GASON)&&defined(GASEXTRA)
        ///store entropy, useful for searching for shocks and related entropy structures
        DoublePos_t entropy;
#endif
        //@}

#ifdef EXTENDEDFOFINFO
	    ///Variables to store extra info
        ///used by VELOCIraptor group finding
        //@{
        ///file id of input file containing particle
    	Int_t oFile;
        ///index of particle in input file
	    Int_t oIndex;
        ///mpi thread id of thread which originally read particle
	    Int_t oTask;
        ///group id (halo, subhalo, etc)
	    Int_t idStruct;
        ///3DFOF envelop id
	    Int_t idFOFHost;
        ///host id
	    Int_t idHost;
        //}
#endif


        public:

        /// \name Constructors & Destructors
        //@{
        Particle(Double_t Mass = 0, Double_t x = 0, Double_t y = 0, Double_t z = 0,
                Double_t vx = 0, Double_t vy = 0, Double_t vz = 0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Int_t PID=0);
        Particle(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Int_t PID=0);
        Particle(const Particle &p);
        Particle(std::istream &F);
        //No dynamic allocation, thus destructor not needed.
        ~Particle(){};
        //@}


        /// \name Overloaded operators
        //@{
        Particle &operator=(const Particle &part);
        bool operator==(const Particle &p) const
        {
            int ival=1;
            ival*=((position[0]==p.position[0])&&(position[1]==p.position[1])&&(position[2]==p.position[2])&&
                (velocity[0]==p.velocity[0])&&(velocity[1]==p.velocity[1])&&(velocity[2]==p.velocity[2])&&
                (id==p.id)&&(type==p.type)&&(rho==p.rho)&&(phi==p.phi)&&
                (pid==p.pid));

#ifndef NOMASS
            ival*=((mass==p.mass));
#endif
#ifdef GASON
            ival*=((u==p.u)&&(sphden==p.sphden));
#endif
#ifdef STARON
            ival*=((tage==p.tage));
#endif
#if defined(GASON) && defined(STARON)
            ival*=((zmet==p.zmet)&&(sfr==p.sfr));
#endif

#if defined(GASON) && (GASEXTRA)
            ival*=((entropy==p.entropy));
#endif
            return ival;
        }
        bool operator!=(const Particle &p) const {return !((*this)==p);}
        // Inserter: prints out the particle in ascii format,
        // e.g., cout << part;
        friend ostream &operator << (ostream &outs, const Particle &p);
        // Extractor: sets particle attributes from a stream,
        // e.g., cin >> part;
        //friend istream &operator >> (istream &ins, Particle &p);
        //@}

        /// \name IO methods
        //@{
        //// Read a particle from an stream in binary format
        //void Read(std::istream &F);
        /// Write a particle to a stream in binary format
        void Write(std::ostream &F) const;
        //@}

        /// \name Get & Set methods
        //@{
#ifdef NOMASS
        Double_t GetMass() const { return MASSVAL; }
        void SetMass(const Double_t &Mass) {}
#else
        Double_t GetMass() const { return mass; }
        void SetMass(const Double_t &Mass) {mass=Mass;}
#endif

#ifndef LOWPRECISIONPOS
        const Double_t* GetPosition() const { return position; }
#else
        const DoublePos_t* GetPosition() const { return position; }
#endif
        Double_t GetPosition(const int &i) const { return position[i]; }
        Double_t X() const { return position[0]; }
        Double_t Y() const { return position[1]; }
        Double_t Z() const { return position[2]; }
        void SetPosition(const int &i, const Double_t &x){ position[i] = x; }
        void SetPosition(const Double_t* p)
        {
        position[0] = p[0];
        position[1] = p[1];
        position[2] = p[2];
        }
        void SetPosition(const Double_t &x, const Double_t &y, const Double_t &z)
        {
        position[0] = x;
        position[1] = y;
        position[2] = z;
        }

#ifndef LOWPRECISIONPOS
        const Double_t *GetVelocity() const { return velocity; }
#else
        const DoublePos_t* GetVelocity() const { return velocity; }
#endif
        Double_t GetVelocity(const int &i) const { return velocity[i]; }
        Double_t Vx() const { return velocity[0]; }
        Double_t Vy() const { return velocity[1]; }
        Double_t Vz() const { return velocity[2]; }
        void SetVelocity(const int &i, const Double_t &x){ velocity[i] = x; }
        void SetVelocity(const Double_t* v)
        {
        velocity[0] = v[0];
        velocity[1] = v[1];
        velocity[2] = v[2];
        }
        void SetVelocity(const Double_t &vx, const Double_t &vy, const Double_t &vz)
        {
        velocity[0] = vx;
        velocity[1] = vy;
        velocity[2] = vz;
        }
        Double_t GetPhase(const int &i){
            if (i<3) return position[i];
            else return velocity[i-3];
        }

        void SetPhase(const int &i, const Double_t &x){
            if (i<3) position[i]=x;
            else velocity[i-3]=x;
        }

#ifdef LONGPIDS
        UInt_t GetPID() {return pid;}
        void SetPID(const UInt_t &i){pid=i;}
#else
        Int_t GetPID() {return pid;}
        void SetPID(const Int_t &i){pid=i;}
#endif
#ifdef LONGIDS
        UInt_t GetID() {return id;}
        void SetID(const UInt_t &i){id=i;}
#else
        Int_t GetID() {return id;}
        void SetID(Int_t i){id=i;}
#endif
        int GetType() {return type;}
        void SetType(int i){type=i;}
        Double_t GetDensity() {return rho;}
        void SetDensity(const Double_t &Rho){rho=Rho;}
        Double_t GetPotential() {return phi;}
        void SetPotential(const Double_t &Phi){phi=Phi;}
#ifdef GASON
        Double_t GetU() {return u;}
        void SetU(const Double_t &U){u=U;}
        Double_t GetSPHDen() {return sphden;}
        void SetSPHDen(const Double_t &SPHden){sphden=SPHden;}
        ///\todo having temperature and entropy functions are not trivial
        ///because need to include eos of gas, metallicity, units, etc.
#endif
#ifdef STARON
        Double_t GetTage() {return tage;}
        void SetTage(const Double_t &Tage){tage=Tage;}
#endif
#if defined (GASON) && (STARON)
        Double_t GetZmet() {return zmet;}
        void SetZmet(const Double_t &Zmet){zmet=Zmet;}
        Double_t GetSFR() {return sfr;}
        void SetSFR(const Double_t &SFR){sfr=SFR;}
#endif

#if defined(GASON) && (GASEXTRA)
        Double_t GetEntropy() {return entropy;}
        void SetEntropy(const Double_t &Entropy) {entropy=Entropy;}
#endif

#ifdef EXTENDEDFOFINFO
        ///Sets and Gets for ExtendedOutput variables
        void SetOFile(const Int_t &i) {oFile = i;}
        Int_t GetOFile() {return oFile;}

        void SetOIndex(const Int_t &i) {oIndex = i;}
        Int_t GetOIndex() {return oIndex;}

        void SetOTask(const Int_t &i) {oTask = i;}
        Int_t GetOTask() {return oTask;}

        void SetIdStruct(const Int_t &i) {idStruct = i;}
        Int_t GetIdStruct() {return idStruct;}

        void SetIdHost(const Int_t &i) {idHost = i;}
        Int_t GetIdHost() {return idHost;}

        void SetIdFOFHost(const Int_t &i) {idFOFHost = i;}
        Int_t GetIdFOFHost() {return idFOFHost;}
#endif

        //@}

        /// \name Other useful functions
        //@{
        ///scale positions by a quantity
        void ScalePositions(Double_t &x){position[0]*=x;position[1]*=x;position[2]*=x;}
        ///scale positions by a quantity
        void ScaleVelocity(Double_t &x){velocity[0]*=x;velocity[1]*=x;velocity[2]*=x;}
        ///scale positions/velocity by a quantity
        void ScalePhase(Double_t &x, Double_t &v){position[0]*=x;position[1]*=x;position[2]*=x;velocity[0]*=v;velocity[1]*=v;velocity[2]*=v;}
        /// Radius of particle
        Double_t Radius() const {return sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);}
        ///Radius squared
        Double_t Radius2() const{return position[0] * position[0]+position[1]*position[1]+position[2]*position[2];}
        /// Spherical coordinate theta (angle between the position vector and the z-axis)
        Double_t Theta() const
        {
            Double_t r=Radius();
            if (r>0) return acos(position[2] / Radius());
            else return 0;
        }
        /// Spherical coordinate phi (azimuthal angle)
        Double_t Phi() const
        {
            if (position[1]>=0&&position[0]>=0)
            return atan2(position[1], position[0]);
            else if (position[1]>0&&position[0]<0)
            return -atan2(position[1], -position[0])+1.0*M_PI;
            else if (position[1]<0&&position[0]<0)
            return atan2(-position[1], -position[0])+1.0*M_PI;
            else
            return -atan2(-position[1], position[0])+2.0*M_PI;
        }
        /// Cylindrical radius (with symmetry axis along z)
        Double_t CylRadius() const
        {
            return sqrt(position[0] * position[0] + position[1] * position[1]);
        }
        /// Absolute velocity of particle
        Double_t AbsoluteVelocity() const
        {
        return sqrt(velocity[0] * velocity[0] +
            velocity[1] * velocity[1] +
            velocity[2] * velocity[2]);
        }
        /// Velocity along radial direction, check if returns nan then
        Double_t RadialVelocity() const
        {
        Double_t r = Radius();
        Double_t rv= (velocity[0] * position[0] / r +
            velocity[1] * position[1] / r +
            velocity[2] * position[2] / r);
        if (std::isnan(rv)) return AbsoluteVelocity();
        else return rv;
        }
        /// Velocity tangential to radial direction check if returns nan
        Double_t CircularVelocity() const
        {
            Double_t r = Radius();
            Double_t x = position[1] * velocity[2] - position[2] * velocity[1];
            Double_t y = position[2] * velocity[0] - position[0] * velocity[2];
            Double_t z = position[0] * velocity[1] - position[1] * velocity[0];
            Double_t cv =sqrt(x * x + y * y + z * z) / r;
            if (std::isnan(cv)) return 0.0;
            else return cv;
        }
        /// Radial velocity in cylindrical coordinates.
        Double_t CylRadialVelocity() const
        {
            Double_t p = Phi();
            return velocity[0] * cos(p) + velocity[1] * sin(p);
        }
        /// Tangential velocity in cylindrical coordinates.
        Double_t CylTangentialVelocity() const
        {
            Double_t p = Phi();
            return -velocity[0] * sin(p) + velocity[1] * cos(p);
        }
        /// Angular momentum, J = r * Vt
        Double_t AngularMomentum() const
        {
            return CircularVelocity() * Radius();
        }
        //@}


#ifdef USEBOOSTMPI
        /// \name Serialization of class needed for mpi using boost library
        //@{
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
#ifndef NOMASS
			ar & mass;
#endif
            ar & position;
            ar & velocity;
            ar & id;
            ar & type;
            ar & rho;
        }
        //@}
#endif

    };

    /*
        End of Particle Class
    */


/*!
    \class NBody::GasParticle
    \brief An sph gas particle

    The class is a subclass of \ref NBody::Particle and has several extra quantities commonly need or used in hydro simulations
    These are temperature, log of entropy, pressure, electron or ionization fraction, atomic hydrogen, metallicity, star foramtion rate, internal energy. \n
    This class is in a state of flux as more properties could be added.

*/
    class GasParticle : public Particle
    {
        private:
        /// \name temperature, log of entropy, pressure, electron or ionization fraction, atomic hydrogen, metallicity, star foramtion rate, internal energy what else
        //@{
        Double_t temp, lgS, P, Ne, Nh0, metal, sfr, U;
        //@}
        public:

        /// \name Constructors & Destructors
        //@{
        GasParticle(Double_t Mass = 0, Double_t x = 0, Double_t y = 0, Double_t z = 0,
                Double_t vx = 0, Double_t vy = 0, Double_t vz = 0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0,
                Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0,Double_t NH0=0,Double_t Zi=0, Double_t SFR=0,Double_t LGS=0);
        GasParticle(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0,
                    Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0,Double_t NH0=0,Double_t Zi=0, Double_t SFR=0,Double_t LGS=0);
        GasParticle(const GasParticle &p);
        //@}

        /// \name Overloaded operators
        //@{
        GasParticle &operator=(const GasParticle &part);
        bool operator==(const GasParticle &p) {return ((Particle::operator==((Particle)p))&&temp==p.temp&&lgS==p.lgS&&P==p.P&&Ne==p.Ne&&Nh0==p.Nh0&&metal==p.metal&&sfr==p.sfr&&U==p.U);}
        bool operator!=(const GasParticle &p) const {return !((Particle::operator==((Particle)p))&&temp==p.temp&&lgS==p.lgS&&P==p.P&&Ne==p.Ne&&Nh0==p.Nh0&&metal==p.metal&&sfr==p.sfr&&U==p.U);}
        //@}

        /// \name Get & Set methods
        //@{
        Double_t GetTemp() const { return temp; }
        void SetTemp(Double_t Temp) { temp = Temp; }
        Double_t GetPressure() const { return P; }
        void SetPressure(Double_t Pi) { P = Pi; }
        Double_t GetU() const { return U; }
        void SetU(Double_t Ui) { U = Ui; }
        Double_t GetNe() const { return Ne; }
        void SetNe(Double_t NE) { Ne = NE; }
        Double_t GetNh0() const { return Nh0; }
        void SetNh0(Double_t NH0) { Nh0 = NH0; }
        Double_t GetZ() const { return metal; }
        void SetZ(Double_t Zi) { metal = Zi; }
        Double_t Getsfr() const { return sfr; }
        void Setsfr(Double_t SFR) { sfr = SFR; }
        Double_t GetEntropy() const { return lgS; }
        void SetEntropy(Double_t lgs) { lgS = lgs; }
        //@}

#ifdef USEBOOSTMPI
        /// \name Serialization of class needed for mpi using boost library
        //@{
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Particle>(*this);
            ar & temp;
            ar & lgS;
            ar & P;
            ar & Ne;
            ar & Nh0;
            ar & metal;
            ar & sfr;
            ar & U;
        }
        //@}
#endif

    };

/*!
    \class NBody::StarParticle
    \brief An sph star particle

    The class is a subclass of \ref NBody::Particle and has several extra quantities commonly need or used in hydro simulations
    These are formation time and metallicity. \n
    This class is in a state of flux as more properties could be added.

*/
    class StarParticle : public Particle
    {
        private:

        ///formation time
        Double_t tform;
        ///metallicity
        Double_t metal;

        public:
        /// \name Constructors & Destructors
        //@{
        StarParticle(Double_t Mass = 0, Double_t x = 0, Double_t y = 0, Double_t z = 0,
                Double_t vx = 0, Double_t vy = 0, Double_t vz = 0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0);
        StarParticle(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0);
        StarParticle(const StarParticle &p);
        //@}

        /// \name Overloaded operators
        //@{
        StarParticle &operator=(const StarParticle &part);
        bool operator==(const StarParticle &p) {return ((Particle::operator==((Particle)p))&&tform==p.tform&&metal==p.metal);}
        bool operator!=(const StarParticle &p) const { return !((Particle::operator==((Particle)p))&&tform==p.tform&&metal==p.metal);}
        //@}

        /// \name Get & Set methods
        //@{
        Double_t GetFormationTime() const { return tform; }
        void SetFormationTime(Double_t Tform) { tform = Tform; }
        Double_t GetZ() const { return metal; }
        void SetZ(Double_t zz) { metal = zz; }
        //@}

#ifdef USEBOOSTMPI
        /// \name Serialization of class needed for mpi using boost library
        //@{
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Particle>(*this);
            ar & tform;
            ar & Z;
        }
        //@}
#endif

    };
    /*
        End of Particle SubClasses
    */

    ///\name useful comparison functions for qsort involving particles
    //@{
    ///sort in ascending particle pid
    int PIDCompare (const void *a, const void *b);
    ///sort in ascending particle id
    int IDCompare (const void *a, const void *b);
    ///sort in ascending particle radius
    int RadCompare (const void *a, const void *b);
    ///sort in ascending particle type
    int TypeCompare (const void *a, const void *b);
    ///sort in ascending particle potential
    int PotCompare (const void *a, const void *b);
    //@}
}

#ifdef USEBOOSTMPI
        namespace boost { namespace mpi {
            template <>
            struct is_mpi_datatype<NBody::Particle> : mpl::true_ { };
        } }
        namespace boost { namespace mpi {
            template <>
            struct is_mpi_datatype<NBody::GasParticle> : mpl::true_ { };
        } }
        namespace boost { namespace mpi {
            template <>
            struct is_mpi_datatype<NBody::StarParticle> : mpl::true_ { };
        } }
#endif

#endif // PARTICLE_H
