/*! \file allvars.h
 *  \brief declares global variables and structures.
 *
 *  This file declares all global variables and structures. Further variables should be added here, and declared as
 *  \e \b extern. The actual existence of these variables is provided by the file \ref allvars.cxx. To produce
 *  \ref allvars.cxx from \ref allvars.h, do the following:
 *
 *     \arg Erase all \#define's, typedef's, and enum's
 *     \arg add \#include "allvars.h", delete the \#ifndef ALLVARS_H conditional
 *     \arg delete all keywords 'extern'
 *     \arg delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/timeb.h>

//gsl stuff
#include <gsl/gsl_heapsort.h>

//nbody code
#include <NBody.h>
#include <NBodyMath.h>
#include <KDTree.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#ifdef USEMPI
#include <mpi.h>
///Includes external global variables used for MPI version of code
#include "mpivar.h"
#endif

#include "../../src/gadgetitems.h"
#include "../../src/tipsy_structs.h"

//for morphological classifications (where tidal tails, or completely disrupted)
#define NUMMORPHCLASS 3;
#define QTIDAL 0.60
#define STIDAL 0.40
#define QLIMTIDAL 0.30
#define SLIMTIDAL 0.25
#define VFTIDAL 0.85
#define VFLIMTIDAL 0.50
#define V2TIDAL 0.6
#define V3TIDAL 0.4
#define V2LIMTIDAL 0.35
#define V3LIMTIDAL 0.30

#define ETIDAL 0.90
#define ELIMTIDAL 0.50

//number of bands in which to calculate magnitude
#define NBANDS 9

///define the particle types
//@{
#define NPARTTYPES 5
#define ALLTYPE 4
#define GASTYPE 0
#define DMTYPE 1
#define STARTYPE 2
#define BHTYPE 3
//@}

///\todo need to alter how stats are stored if one can update the number of quantiles stored
///the number of values that is used to characterize a bulk distribution
///where idea is 0=mean,1=sd,2-4=quantiles of 0.16,0.5,0.84
#define NDISTRIB 5
///max number of quantiles
#define MAXNQUANTS 11

///max number of radial bins
#define NRAD 1000


/// Structures and external variables

///for stars,gas, define lower limits for Z and sfr
#define MINGASZ 1e-7
#define MINSTARZ 1e-7
#define MINSFR 1e-3

///define sphi io block points for gadget
//@{
///self energy  
#define sphioblock_u 0
///density
#define sphioblock_d 1
#ifdef SPHCOOLING
///mean molecular weight
#define sphioblock_mmw 2
///neutral hydrogen abundance
#define sphioblock_nha 3
#endif
///smoothing lengths
#define sphioblock_h
//@}

///define different reference frames
//@{
///cm ref
#define ICMFRAME 1
///particle with largest potential
#define IPOTFRAME 2
///most bound particle based in input data
#define IMBFRAME 3
//@}

///define omp number of elements threshold
#define ompunbindnum 1000
#define omppropnum 50000

///minimum number of particles for properties to be calculated
#define MINNPROP 10 
///minimum number of particles for a profile to be calculated
#define MINNPROFILE 100
///minimum mass fraction of candidate structure for enclosed quantities to be useful
#define MINMASSFRAC 0.05

///some useful constants for gas
//@{
///mass of helium relative to hydrogen
#define M_HetoM_H 4.0026
//@}

using namespace std;
using namespace Math;
using namespace NBody;

///for potential calculation
struct PotInfo
{
    ///bucket size for tree calculation
    int BucketSize;
    ///for tree potential calculation;
    Double_t TreeThetaOpen;
    ///gravitational softening (here simple plummer)
    Double_t eps;
    ///\todo When calculating the thermal internal energy of gas particles
    ///I should probably include some temperature floors and ceilings
    ///and perhaps some sound speed constrains. But I should produce
    ///properties of the baryonic material for all, gravitationally self-bound and grav+therm bound

    PotInfo(){
        BucketSize=8;
        TreeThetaOpen=0.7;
        eps=0.01;
    }
};

///Options structure
struct Options
{
    ///filenames
    //@{
    ///input particle data and output name
    char *fname,*outname;
    ///halo catalog base name
    char *halocatname;
    ///optional config file name
    char *configname;
    ///optional stellar synthesis model file name
    char *magname;
    //@}
    ///number of files per particle data snapshot, if zero then tipsy file
    int nfiles;
    ///number of files per halo catalog
    int nfileshalocat;
    ///number of sph blocks in gadget file
    int sphnio;
    ///defining units and cosmology
    //@{
    ///length,m,v,grav units pressure unit, temperature
    Double_t L, M, V, G, PUnit, TUnit;
    ///scalbodydata/utils/haloanalcode/e factor, Hubunit, h and cosmology
    Double_t a,H,h, Omega_m, Omega_Lambda,rhoc;
    ///softening length and period (comove)
    Double_t eps, p;
    ///background density and virial level
    Double_t rhobg, virlevel;
    int comove;
    ///for gas, mass fraction of hydrogen
    Double_t Xh;
    ///for gas, boltzmann constant
    Double_t Boltzmann;
    ///for gas, number of hydrogen atoms per sph particle
    Double_t N_H;
    ///for gas, correction for helium to mass of sph particle /mh ie: Number of particles = M_p / (M_H*Xh+M_He*(1-Xh)) or even M_H*Xh+M_He*(1-Xh)*Y_Li+M_Li*(1-Xh)*(1-Y_Li) ...
    Double_t N_He_corr;
    ///for gas, relation of specific internal energy to Temperature (u=3/2*k_B*(N_H*(Xh+(1-Xh)*N_He)*T
    Double_t utoT;
    ///for gas, temperature threshold when calculating X-ray based spectroscopic-like temperature)
    Double_t Tempthreshold;
    ///\todo may want some specific stuff for stars too
    //@}

    ///store number of particles
    Int_t npart[NPARTTYPES];

    ///structure that contains variables for potentialcalc
    PotInfo pinfo;
    Double_t error;

    ///quantile points for characterization of distribution
    int nquants;
    Double_t quants[MAXNQUANTS];

    ///minimum number a halo has to have to be analyzed
    Int_t minnum;
    ///Kinetic/Potential Energy ratio
    Double_t TVratio;
    ///number of mpi files that were written in parallel that must be read for group_catalog files
    int mpinum;

    ///for calculaion of CM
    Double_t cmfrac, cmadjustfac;
    PotInfo uinfo;

    ///for calculating kinetic energy to a reference frame
    int reframe;

    ///if no mass is stored
    Double_t MassValue;

    ///if zoom simulations, store effective resolution
    long int Neff;

    ///for verbose output
    int verbose;

    ///store whether potential has been already calculated as object has been moved to particle with deepest potential well. 
    int ipotcalc;

    ///used to calculate kinetic reference frame about particle at bottom of potential well
    Int_t Nbref;
    Double_t fracbref;

    ///input STF/VELOCIraptor format flags
    //@{
    int iseparatefiles,ibinaryin;
    //@}

    ///output format
    int ibinaryout;

    Options()
    {
        fname=outname=halocatname=NULL;
        nfiles=1;
        sphnio=2;
        nfileshalocat=1;

        L = 1.0;
        M = 1.0;
        V = 1.0;
        PUnit=5.588e-14;//(assumes T in kelvin) makes pressure into ergs/cm^3 assuming standard gadget units (this is boltzmann's constant / proton mass *Massunit_g/(Lengthunit_cm)^3
        TUnit=1.0;
        MassValue=1.0;

        p = 0.0;
        a = 1.0;
        H = 0.1;
        h = 1.0;
        Omega_m = 1.0;
        Omega_Lambda = 0.0;
        rhobg = 1.0;
        virlevel = 50.;
        comove=0;
        G = 43.01;//for G in units of (km/s)^2 Mpc/1e10 Msun
        rhoc=27.753;//for same units as above, of course missing factors of h^2
        G = 43021.1349;//for G in units of (km/s)^2 kpc/1e10 Msun
        rhoc=2.7753e-8;

        Xh=0.76;
        Boltzmann=log(6.94107285)-70.0*log(10.0); //ln of boltzmann in units of K^(-1) (km/s)^2 (1e10 Msun)
        N_H=log(1.1892)+67.0*log(10.0); //number of hydrogen atoms per 10^10 Msun
        N_He_corr=-log(Xh+M_HetoM_H*(1.0-Xh));
        utoT=log(1.5)+Boltzmann+N_H+N_He_corr; 
        Tempthreshold=7.0*log(10.0);

        error=1e-3;

        nquants=3;
        quants[0]=0.16;quants[1]=0.5;quants[2]=0.84;
        for (int i=0;i<NPARTTYPES;i++) npart[i]=0;

        minnum=20;
        TVratio=-1.0;
        mpinum=0;

        cmfrac=0.1;
        cmadjustfac=0.7;
        Neff=-1;

        reframe=IPOTFRAME;
        verbose=0;

        ipotcalc=1;
        Nbref=10;
        fracbref=0.05;

        iseparatefiles=0;
        ibinaryin=0;
        ibinaryout=0;

    }
};

///Radial profile data in spherical shells
struct RadPropData
{
    ///To store radial profiles
    //@{
    Int_t nbins;
    int ptype;
    Double_t *numinbin;
    Double_t *radval;
    Double_t *Menc;
    Coordinate *velenc;
    Matrix *sigmavenc;
    Coordinate *Jenc;
    Double_t *qenc;
    Double_t *senc;
    Matrix *eigvecenc;
    ///primarily for sph particles
    Double_t *Den;
    ///gas/star particle info
    Double_t *Tempenc;
    Double_t *Tdispenc;
    Double_t *Temenc;
    Double_t *Tslenc;
#ifdef RADIATIVE
    Double_t *Zmetenc;
    Double_t *Tageenc;
#endif
//@}
    RadPropData(){
        //for (int i=0;i<NRAD;i++) {numinbin[i]=radval[i]=Menc[i]=0.;Jenc[i]=Coordinate(0,0,0);velenc[i]=Coordinate(0,0,0);qenc[i]=senc[i]=1.0;}
        ptype=DMTYPE;
        numinbin=NULL;radval=NULL;Menc=NULL;velenc=NULL;sigmavenc=NULL;Jenc=NULL;qenc=NULL;senc=NULL;eigvecenc=NULL;
        Den=NULL;Tempenc=NULL;Tdispenc=NULL;Temenc=NULL;Tslenc=NULL;
#ifdef RADIATIVE
        Zmetenc=NULL;
        Tageenc=NULL;
#endif
    }
    void SetBins(int N) {
        nbins=N;
        if(numinbin!=NULL) delete[] numinbin;numinbin=new Double_t[nbins];
        if(radval!=NULL) delete[] radval;radval=new Double_t[nbins];
        if(Menc!=NULL) delete[] Menc;Menc=new Double_t[nbins];
        if(velenc!=NULL) delete[] velenc;velenc=new Coordinate[nbins];
        if(sigmavenc!=NULL) delete[] sigmavenc;sigmavenc=new Matrix[nbins];
        if(Jenc!=NULL) delete[] Jenc;Jenc=new Coordinate[nbins];
        if(qenc!=NULL) delete[] qenc;qenc=new Double_t[nbins];
        if(senc!=NULL) delete[] senc;senc=new Double_t[nbins];
        if(eigvecenc!=NULL) delete[] eigvecenc;eigvecenc=new Matrix[nbins];
        if(Tempenc!=NULL) delete[] Tempenc;Tempenc=new Double_t[nbins];
        if(Tdispenc!=NULL) delete[] Tdispenc;Tdispenc=new Double_t[nbins];
        if(Temenc!=NULL) delete[] Temenc;Temenc=new Double_t[nbins];
        if(Tslenc!=NULL) delete[] Tslenc;Tslenc=new Double_t[nbins];
        if(Den!=NULL) delete[] Den;Den=new Double_t[nbins];
#ifdef RADIATIVE
        if(Zmetenc!=NULL) delete[] Zmetenc;Zmetenc=new Double_t[nbins];
        if(Tageenc!=NULL) delete[] Tageenc;Tageenc=new Double_t[nbins];
#endif
        for (int i=0;i<nbins;i++) {numinbin[i]=radval[i]=Menc[i]=0;Jenc[i]=Coordinate(0,0,0);velenc[i]=Coordinate(0,0,0);qenc[i]=senc[i]=1.0;}
        for (int i=0;i<nbins;i++) {Tempenc[i]=Tdispenc[i]=Den[i]=Temenc[i]=Tslenc[i]=0.;}
    }
};

///Cylindrical profile data
struct CylPropData
{
    ///reference axis of cylinder
    Coordinate zref;
    ///To store cylindrical profiles 
    //@{
    Int_t nbins;
    Double_t *numinbin;
    Double_t *radval;
    Double_t *Menc;
    Coordinate *Jenc;
    Coordinate *velenc;
    Matrix *sigmavenc;
    Double_t *zmean;
    Double_t *Rmean;
    Double_t *zdisp;
    Double_t *Rdisp;
    Double_t *qenc;
    Double_t *senc;
    Matrix *eigvecenc;
    ///primarily for sph particles
    Double_t *Den;
    ///gas/star particle info
    Double_t *Tempenc;
    Double_t *Tdispenc;
    Double_t *Temenc;
    Double_t *Tslenc;
#ifdef RADIATIVE
    Double_t *Zmetenc;
    Double_t *Tageenc;
#endif
    //@}
    CylPropData(){
        numinbin=NULL;radval=NULL;Menc=NULL;velenc=NULL;sigmavenc=NULL;Jenc=NULL;zmean=NULL;Rmean=NULL;zdisp=NULL;Rdisp=NULL;qenc=NULL;senc=NULL;eigvecenc=NULL;
        Den=NULL;Tempenc=NULL;Tdispenc=NULL;Temenc=NULL;Tslenc=NULL;
#ifdef RADIATIVE
        Zmetenc=NULL;
        Tageenc=NULL;
#endif
    }
    void SetBins(int N) {
        nbins=N;
        if(numinbin!=NULL) delete[] numinbin;numinbin=new Double_t[nbins];
        if(radval!=NULL) delete[] radval;radval=new Double_t[nbins];
        if(Menc!=NULL) delete[] Menc;Menc=new Double_t[nbins];
        if(velenc!=NULL) delete[] velenc;velenc=new Coordinate[nbins];
        if(sigmavenc!=NULL) delete[] sigmavenc;sigmavenc=new Matrix[nbins];
        if(Jenc!=NULL) delete[] Jenc;Jenc=new Coordinate[nbins];
        if(zmean!=NULL) delete[] zmean;zmean=new Double_t[nbins];
        if(Rmean!=NULL) delete[] Rmean;Rmean=new Double_t[nbins];
        if(zdisp!=NULL) delete[] zdisp;zdisp=new Double_t[nbins];
        if(Rdisp!=NULL) delete[] Rdisp;Rdisp=new Double_t[nbins];
        if(qenc!=NULL) delete[] qenc;qenc=new Double_t[nbins];
        if(senc!=NULL) delete[] senc;senc=new Double_t[nbins];
        if(eigvecenc!=NULL) delete[] eigvecenc;eigvecenc=new Matrix[nbins];
        if(Tempenc!=NULL) delete[] Tempenc;Tempenc=new Double_t[nbins];
        if(Tdispenc!=NULL) delete[] Tdispenc;Tdispenc=new Double_t[nbins];
        if(Temenc!=NULL) delete[] Temenc;Temenc=new Double_t[nbins];
        if(Tslenc!=NULL) delete[] Tslenc;Tslenc=new Double_t[nbins];
        if(Den!=NULL) delete[] Den;Den=new Double_t[nbins];
#ifdef RADIATIVE
        if(Zmetenc!=NULL) delete[] Zmetenc;Zmetenc=new Double_t[nbins];
        if(Tageenc!=NULL) delete[] Tageenc;Tageenc=new Double_t[nbins];
#endif
        for (int i=0;i<nbins;i++) {numinbin[i]=radval[i]=Menc[i]=0.;Jenc[i]=velenc[i]=Coordinate(0,0,0);qenc[i]=senc[i]=1.0;zmean[i]=Rmean[i]=zdisp[i]=Rdisp[i]=0.;}
        for (int i=0;i<nbins;i++) {Tempenc[i]=Tdispenc[i]=Den[i]=Temenc[i]=Tslenc[i]=0.;}
    }
};

///Based on PropData in \ref allvars.h of STF/VELOCIraptor
///Here the data also stores the type of particles it corresponds to
struct PropData
{
    ///type of structure, dm, gas, star, etc
    Int_t ptype;
    ///number of particles
    Int_t num;
    ///centre of mass
    Coordinate gcm, gcmvel;
    ///Position of most bound particle
    Coordinate gpos, gvel;
    ///physical properties regarding mass, size
    Double_t gmass,gsize,gMvir,gRvir,gRmbp,gmaxvel,gRmaxvel,gMmaxvel;
    ///physical properties for shape/mass distribution
    Double_t gq,gs;
    Matrix geigvec;
    ///physical properties for velocity
    Double_t gsigma_v;
    Matrix gveldisp;
    ///physical properties for dynamical state
    Double_t Efrac,Pot,T;
    ///physical properties for dynamical state for specific particle types
    Double_t Efractyped[NPARTTYPES],Pottyped[NPARTTYPES],Ttyped[NPARTTYPES];
    ///physical properties for angular momentum
    Coordinate gJ;

    ///metallicity, temperature, pressure,etc
    Double_t Temp[NDISTRIB];
#ifdef RADIATIVE
    Double_t Zmet[NDISTRIB];
    ///stellar age
    Double_t Tage[NDISTRIB];
    ///luminosity, which for gas which would probably be Infrared and radio bands
    Double_t LUM[NBANDS];
#endif

    ///To store cylindrical profiles 
    RadPropData radprofile;
    ///To store cylindrical profiles 
    CylPropData cylprofile;


    PropData(){
        num=0;
        gmass=gsize=gRmaxvel=Efrac=Pot=T=gcm[0]=gcm[1]=gcm[2]=gcmvel[0]=gcmvel[1]=gcmvel[2]=0.;
        gveldisp=Matrix(0.);
        gq=gs=1.0;
        gRmbp=0.0;

        for (int i=0;i<NDISTRIB;i++)Temp[i]=0;
#ifdef RADIATIVE
        for (int i=0;i<NDISTRIB;i++)Zmet[i]=Tage[i]=0.;
        for (int i=0;i<NBANDS;i++)LUM[i]=0;
#endif
    }
};


///used to store halo particle id data
struct HaloParticleData{
    ///store the halo id
    long unsigned haloID;
    ///number of particles, eventually used to store bound number
    long unsigned NumberofParticles;
    ///number of particles of different types, eventually used to store bound number
    Int_t NumofType[NPARTTYPES];
    ///number of all particles
    long unsigned AllNumberofParticles;
    ///number of all particles of different types
    Int_t AllNumofType[NPARTTYPES];

    ///offset of halo's particle in list of all particles in structures
    Int_t noffset;
    ///array of particle ids
    long unsigned *ParticleID;
    HaloParticleData(long unsigned numinhalo=0){
        AllNumberofParticles=NumberofParticles=numinhalo;
        if (NumberofParticles>0) {
            ParticleID=new long unsigned[NumberofParticles];
        }
        for (int i=0;i<NPARTTYPES;i++) AllNumofType[i]=NumofType[i]=0;
    }
    void Alloc(long unsigned ninhalos=0){
        if (AllNumberofParticles>0){
            delete[] ParticleID;
        }
        AllNumberofParticles=NumberofParticles=ninhalos;
        if (AllNumberofParticles>0) {
            ParticleID=new long unsigned[NumberofParticles];
        }
    }
    ~HaloParticleData(){
        if (NumberofParticles>0){
            delete[] ParticleID;
        }
    }
};

///External types for tipsy
//@{
#define TGASTYPE 0;
#define TDARKTYPE 1;
#define TSTARTYPE 2;
//@}

///note that for grouped particles type is 10+TYPE.
#define SUBSTRUCTTYPE 10;

#endif

