/*! \file allvars.h
 *  \brief declares global variables.
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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <getopt.h>
#include <sys/stat.h> 
#include <sys/timeb.h>

#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


///\name Include for NBodyFramework library.
//@{
///nbody code
#include <NBody.h>
///Math code
#include <NBodyMath.h>
///Binary KD-Tree code
#include <KDTree.h>
///Extra routines that analyze a distribution of particles
#include <Analysis.h>
//@}

// ///\name for checking the endian of floats
//#include <endianutils.h>

///if using OpenMP API
#ifdef USEOPENMP
#include <omp.h>
#endif

///if using HDF API
#ifdef USEHDF
#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#endif

using namespace std;
using namespace Math;
using namespace NBody;

//-- Structures and external variables

/// \name Particle types
//@{
/// \defgroup PARTTYPES
//@{
#define  GASTYPE 0
#define  DARKTYPE 1
#define  DARK2TYPE 2
#define  DARK3TYPE 3
#define  STARTYPE 4
#define  BHTYPE 5
#define  WINDTYPE 6
#define  NPARTTYPES 7
//number of baryon types +1, to store all baryons
#define  NBARYONTYPES 5
//@}
//@}

/// \defgroup SEARCHTYPES 
//@{
/// \name Specify particle type to be searched, all, dm only, separate
//@{
#define  PSTALL 1 
#define  PSTDARK 2 
#define  PSTSTAR 3 
#define  PSTGAS 4 
#define  PSTBH 5 
//@}
//@}

/// \defgroup STRUCTURETYPES
//@{
/// \name Specific structure type, allow for other types beside HALO
/// \todo note that here I have set background group type to a halo structure type but that can be changed
//@{
#define HALOSTYPE 10
#define HALOCORESTYPE 5
#define WALLSTYPE 1
#define VOIDSTYPE 2
#define FILAMENTSTYPE 3
#define BGTYPE 10
#define GROUPNOPARENT -1
//@}
//@}

/// \defgroup FOFTYPES
//@{
/// \name FOF search types
//@{
//subsets made
///call \ref FOFStreamwithprob
#define  FOFSTPROB 1 
///6D FOF search
#define  FOF6DSUBSET 7
///like \ref FOFStreamwithprob search but search is limited to nearest physical neighbours
#define  FOFSTPROBNN 9
///like \ref FOFStreamwithprob search but here linking length adjusted by velocity offset, smaller lengths for larger velocity offsets
#define  FOFSTPROBLX 10
///like \ref FOFSTPROBLX but for NN search
#define  FOFSTPROBNNLX 11
///like \ref FOFSTPROBNN but there is not linking length applied just use nearest neighbours
#define  FOFSTPROBNNNODIST 12
//for iterative method with FOFStreamwithprob
//#define  FOFSTPROBIT 13
#define  FOFSTPROBSCALEELL 13
//#define  FOFSTPROBIT 13
#define  FOFSTPROBSCALEELLNN 14

//no subsets made, just 6d , 3d and FOFStream searches
#define  FOF6D 4
#define  FOF3D 5
#define  FOFSTNOSUBSET 6
//@}
//@}

/// \defgroup INTERATIVESEARCHPARAMS
//@{
/// \name for iterative subsubstructure search
//@{
/// this is minimum particle number size for a subsearch to proceed whereby substructure split up into CELLSPLITNUM new cells 
#define  MINCELLSIZE 200
#define  CELLSPLITNUM 16
#define  MINSUBSIZE MINCELLSIZE*CELLSPLITNUM
#define  MAXSUBLEVEL 8
//@}
//@}

///\defgroup GRIDTYPES
//@{
/// \name Type of Grid structures 
//@{
#define  PHYSENGRID 1
#define  PHASEENGRID 2
#define  PHYSGRID 3
//@}
//@}

/// \name Max number of neighbouring cells used in interpolation of background velocity field.
//@{

//if cells were cubes this would be all the neighbours that full enclose the cube, 6 faces+20 diagonals
//using daganoals may not be ideal. Furthermore, code using adaptive grid that if effectively produces
//cell that are rectangular prisms. Furthermore, not all cells will share a boundary with another cell.
//So just consider "faces".
#define MAXNGRID 6

//@}

///\defgroup INPUTTYPES
//@{
/// \name defining types of input
//@{
#define  IOGADGET 1
#define  IOHDF 2
#define  IOTIPSY 3
#define  IORAMSES 4
//@}
//@}

/// \name For Unbinding
//@{

///number below which just use PP calculation for potential, which occurs roughly at when n~2*log(n) (from scaling of n^2 vs n ln(n) for PP vs tree and factor of 2 is
///for extra overhead in producing tree. For reasonable values of n (>100) this occurs at ~100. Here to account for extra memory need for tree, we use n=3*log(n) or 150
#define UNBINDNUM 150
///when unbinding check to see if system is bound and least bound particle is also bound
#define USYSANDPART 0
///when unbinding check to see if least bound particle is also bound
#define UPART 1
///use the bulk centre of mass velocity to define velocity reference frame when determining if particle bound
#define CMVELREF 0
///use the particle at potential minimum. Issues if too few particles used as particles will move in and out of deepest point of the potential well
#define POTREF 1

//@}

/// \name For Tree potential calculation
//@{

///leafflag indicating in tree-calculation of potential, reached a leaf node that does not satisfy mono-pole approx
#define leafflag 1
///split flag means node not used as subcells are searched
#define splitflag -1
///cellflag means a node that is not necessarily a leaf node can be approximated by mono-pole
#define cellflag 0

//@}

/// \defgroup OMPLIMS 
//@{
///\name For determining whether loop contains enough for openm to be worthwhile.
//@{
#define ompsearchnum 50000
#define ompunbindnum 1000
#define ompperiodnum 50000
#define omppropnum 50000
//@}
//@}


///\name halo id modifers used with current snapshot value to make temporally unique halo identifiers
#ifdef LONGINT
#define HALOIDSNVAL 1000000000000L
#else
#define HALOIDSNVAL 1000000
#endif

///\defgroup GASPARAMS
//@{
///\name Useful constants for gas
//@{
///mass of helium relative to hydrogen
#define M_HetoM_H 4.0026
//@}
//@}


/// Structure stores unbinding information
struct UnbindInfo 
{
    ///\name flag whether unbind groups, keep bg potential when unbinding, type of unbinding and reference frame
    //@{
    int unbindflag,bgpot,unbindtype,cmvelreftype;
    //@}
    ///fraction of potential energy that kinetic energy is allowed to be and consider particle bound
    Double_t Eratio;
    ///minimum bound mass fraction 
    Double_t minEfrac;
    ///when to recalculate kinetic energies if cmvel has changed enough
    Double_t cmdelta;
    ///maximum fraction of particles to remove when unbinding in one given unbinding step
    Double_t maxunbindfrac;
    ///minimum number of particles to use to calculate reference frame if using particles around deepest potential well as reference frame
    Int_t Npotref;
    ///fraction of number of particles to use to calculate reference frame if using particles around deepest potential well as reference frame
    Double_t fracpotref;
    ///\name gravity and tree potential calculation;
    //@{
    int BucketSize;
    Double_t TreeThetaOpen;
    ///softening length
    Double_t eps;
    //@}
    UnbindInfo(){
        unbindflag=0;
        bgpot=1;
        unbindtype=UPART;
        cmvelreftype=CMVELREF;
        cmdelta=0.02;
        Eratio=1.0;
        minEfrac=1.0;
        BucketSize=8;
        TreeThetaOpen=0.7;
        eps=0.0;
        maxunbindfrac=0.05;
        Npotref=10;
        fracpotref=0.1;
    }
};

/// Structure stores information used when calculating bulk (sub)structure properties
/// which is used in \ref substructureproperties.cxx
struct PropInfo
{
    //interate till this much mass in contained in a spherical region to calculate cm quantities
    Double_t cmfrac,cmadjustfac;

    PropInfo(){
        cmfrac=0.1;
        cmadjustfac=0.7;
    }
};

/// Options structure stores useful variables that have user determined values which are altered by \ref GetArgs in \ref ui.cxx
struct Options
{
    ///\name filenames
    //@{
    char *fname,*outname,*smname,*pname,*gname;
    //@}
    ///input format
    int inputtype;
    ///number of snapshots
    int num_files,snum;
    ///if parallel reading, number of files read in parallel
    int nsnapread;
    ///for output, specify the formats, ie. many separate files, binary or ascii
    int iseparatefiles,ibinaryout;
    ///return propery data in in comoving little h units instead of standard physical units
    int icomoveunit;

    ///\name length,m,v,grav units
    //@{
    Double_t L, M, V, G;
    //@}
    ///period (comove)
    Double_t p;
    ///\name scale factor, Hubunit, h, cosmology, virial density. These are used if linking lengths are scaled or trying to define virlevel using the cosmology
    //@{
    Double_t a,H,h, Omega_m, Omega_b, Omega_cdm, Omega_Lambda, w_de, rhobg, virlevel;
    int comove;
    //@}

    ///to store number of each particle types so that if only searching one particle type, assumes order of gas, dark (halo, disk,bulge), star, special or sink for pfof tipsy style output
    Int_t numpart[NPARTTYPES];

    ///\name parameters that control the local and average volumes used to calculate the local velocity density and the mean field, also the size of the leafnode in the kd-tree used when searching the tree for fof neighbours
    //@{
    int Nvel, Nsearch, Bsize;
    Int_t Ncell;
    Double_t Ncellfac;
    //@}
    ///minimum group size
    int MinSize;
    ///allows for field halos to have a different minimum size
    int HaloMinSize;
    ///Significance parameter for groups
    Double_t siglevel;

    ///whether to search for substructures at all
    int iSubSearch;
    ///type of search
    int foftype,fofbgtype;
    ///grid type, physical, physical+entropy splitting criterion, phase+entropy splitting criterion. Note that this parameter should not be changed from the default value
    int gridtype;
    ///flag indicating search all particle types or just dark matter
    int partsearchtype;
    ///flag indicating a separate baryonic search is run, looking for all particles that are associated in phase-space
    ///with dark matter particles that belong to a structure
    int iBaryonSearch;
    ///flag indicating if move to CM frame for substructure search
    int icmrefadjust;

    ///threshold on particle ELL value, normalized logarithmic distance from predicted maxwellian velocity density.
    Double_t ellthreshold;
    ///\name fofstream search parameters
    //@{
    Double_t thetaopen,Vratio,ellphys;
    //@}
    ///fof6d search parameters
    Double_t ellvel;
    ///scaling for ellphs and ellvel
    Double_t ellxscale,ellvscale;
    ///flag to use iterative method
    int iiterflag;
    ///\name factors used to multiply the input values to find initial candidate particles and for mergering groups in interative search
    //@{
    Double_t ellfac,ellxfac,vfac,thetafac,nminfac;
    Double_t fmerge;
    //@}
    ///factors to alter halo linking length search (related to substructure search)
    //@{
    Double_t ellhalophysfac,ellhalovelfac;
    //@}
    //@{
    ///\name factors used to check for halo mergers, large background substructures and store the velocity scale when searching for associated baryon substructures
    //@{
    Double_t HaloMergerSize,HaloMergerRatio,HaloSigmaV,HaloVelDispScale;
    Double_t fmergebg;
    //@}
    ///flag indicating a single halo is passed or must run search for FOF haloes
    Int_t iSingleHalo;
    ///flag indicating haloes are to be check for self-boundness after being searched for substructure
    Int_t iBoundHalos;
    /// store denv ratio statistics
    //@{
    int idenvflag;
    Double_t denvstat[3];
    //@}
    ///verbose output flag
    int iverbose;
    ///whether or not to write a fof.grp tipsy like array file
    int iwritefof;
    ///whether mass properties for field objects are inclusive
    int iInclusiveHalo;

    ///if no mass value stored then store global mass value
    Double_t MassValue;

    ///structure that contains variables for unbinding
    UnbindInfo uinfo;
    ///structure that contains variables for property calculation
    PropInfo pinfo;

    ///effective resolution for zoom simulations
    Int_t Neff;

    ///\name extra stuff for halo merger check and identification of multiple halo core
    //@{
    int iHaloCoreSearch;
    Double_t halocorexfac, halocorevfac, halocorenfac;
    //@}
    ///for storing a snapshot value to make halo ids unique across snapshots
    long long snapshotvalue;

    ///\name for reading gadget info with lots of extra sph, star and bh blocks 
    //@{
    int gnsphblocks,gnstarblocks,gnbhblocks;
    //@}
    

    ///\name extra runtime flags
    //@{
    ///scale lengths. Useful if searching single halo system and which to automatically scale linking lengths
    int iScaleLengths;
    //@}
    Options()
    {
        L = 1.0;
        M = 1.0;
        V = 1.0;
        G = 1.0;
        p = 0.0;

        a = 1.0;
        H = 0.;
        h = 1.0;
        Omega_m = 1.0;
        Omega_Lambda = 0.0;
        Omega_b = 0.0;
        Omega_cdm = Omega_m;
        rhobg = 1.0;
        virlevel = -1;
        comove=0;
        H=100.0;//value of Hubble flow in h 1 km/s/Mpc
        MassValue=1.0;

        inputtype=IOGADGET;

        num_files=1;
        nsnapread=1;

        fname=outname=smname=pname=gname=outname=NULL;

        Bsize=16;
        Nvel=32;
        Nsearch=256;
        Ncellfac=0.01;

        iSubSearch=1;
        partsearchtype=PSTALL;
        for (int i=0;i<NPARTTYPES;i++)numpart[i]=0;
        foftype=FOFSTPROB;
        gridtype=PHYSENGRID;
        fofbgtype=FOF6D;
        idenvflag=0;
        iBaryonSearch=0;
        icmrefadjust=1;

        Neff=-1;

        ellthreshold=1.5;
        thetaopen=0.05;
        Vratio=1.25;
        ellphys=0.2;
        MinSize=20;
        HaloMinSize=-1;
        siglevel=2.0;
        ellvel=0.5;
        ellxscale=ellvscale=1.0;
        ellhalophysfac=ellhalovelfac=1.0;

        iiterflag=0;
        ellfac=2.5;
        ellxfac=3.0;
        vfac=1.0;
        thetafac=1.0;
        nminfac=0.5;
        fmerge=0.25;

        HaloMergerSize=10000;
        HaloMergerRatio=0.2;
        HaloVelDispScale=0;
        fmergebg=0.5;
        iSingleHalo=0;
        iBoundHalos=0;
        iInclusiveHalo=0;

        iHaloCoreSearch=0;
        halocorexfac=0.5;
        halocorevfac=2.0;
        halocorenfac=0.1;

        iverbose=0;
        iwritefof=0;
        iseparatefiles=0;
        ibinaryout=0;
        icomoveunit=0;

        snapshotvalue=0;

        gnsphblocks=4;
        gnstarblocks=2;
        gnbhblocks=2;

        iScaleLengths=0;
    }
};

/// N-dim grid cell
struct GridCell
{
    int ndim;
    Int_t gid;
    //middle of cell, and boundaries of cell
    Double_t xm[6], xbl[6],xbu[6];
    //mass, radial size of in cell
    Double_t mass, rsize;
    //number of particles in cell
    Int_t nparts,*nindex;
    //neighbouring grid cells and distance from cell centers
    Int_t nnidcells[MAXNGRID];
    Double_t nndist[MAXNGRID];
    Double_t den;
    GridCell(int N=3){
        ndim=N;
        nparts=0;
        den=0;
    }
    ~GridCell(){
        if (nparts>0)delete nindex;
    }
};

/*! structure stores bulk properties like
    \f$ m,\ (x,y,z)_{\rm cm},\ (vx,vy,vz)_{\rm cm},\ V_{\rm max},\ R_{\rm max}, \f$
    which is calculated in \ref substructureproperties.cxx
*/
struct PropData
{
    ///\name order in structure hierarchy and number of subhaloes
    //@{
    long long haloid,hostid,directhostid;
    Int_t numsubs;
    //@}

    ///\name properties of total object including DM, gas, stars, bh, etc
    //@{
    ///number of particles
    Int_t num;
    ///number of particles in FOF envelop
    Int_t gNFOF;
    ///centre of mass
    Coordinate gcm, gcmvel;
    ///Position of most bound particle
    Coordinate gpos, gvel;
    ///\name physical properties regarding mass, size
    //@{
    Double_t gmass,gsize,gMvir,gRvir,gRmbp,gmaxvel,gRmaxvel,gMmaxvel,gRhalfmass;
    Double_t gM200c,gR200c,gM200m,gR200m,gMFOF;
    //@}
    ///\name physical properties for shape/mass distribution
    //@{
    Double_t gq,gs;
    Matrix geigvec;
    //@}
    ///\name physical properties for velocity
    //@{
    Double_t gsigma_v;
    Matrix gveldisp;
    //@}
    ///physical properties for dynamical state
    Double_t Efrac,Pot,T;
    ///physical properties for angular momentum
    Coordinate gJ;
    ///Keep track of position of least unbound particle and most bound particle pid
    UInt_t iunbound,ibound;
    ///Type of structure
    int stype;
    ///concentration (and related quantity used to calculate a concentration)
    Double_t cNFW, VmaxVvir2;
    ///Bullock & Peebles spin parameters
    Double_t glambda_B,glambda_P;
    ///measure of rotational support
    Double_t Krot;
    //@}
#ifdef GASON
    ///\name gas specific quantities
    //@{
    ///number of particles
    int n_gas;
    ///mass
    Double_t M_gas;
    ///pos/vel info
    Coordinate cm_gas,cmvel_gas;
    ///velocity/angular momentum info
    Double_t Krot_gas;
    Coordinate L_gas;
    Matrix veldisp_gas;
    ///morphology
    Double_t Rhalfmass_gas,q_gas,s_gas;
    Matrix eigvec_gas;
    ///mean temperature,metallicty,star formation rate
    Double_t Temp_gas,Z_gas,SFR_gas;
    ///physical properties for dynamical state
    Double_t Efrac_gas,Pot_gas,T_gas;
    //@}
#endif
    
#ifdef STARON
    ///\name star specific quantities
    //@{
    ///number of particles
    int n_star;
    ///mass
    Double_t M_star;
    ///pos/vel info
    Coordinate cm_star,cmvel_star;
    ///velocity/angular momentum info
    Double_t Krot_star;
    Coordinate L_star;
    Matrix veldisp_star;
    ///morphology
    Double_t Rhalfmass_star,q_star,s_star;
    Matrix eigvec_star;
    ///mean age,metallicty
    Double_t t_star,Z_star;
    ///physical properties for dynamical state
    Double_t Efrac_star,Pot_star,T_star;
    //@}
#endif

#ifdef BHON
    ///\name black hole specific quantities
    //@{
    ///number of BH
    int n_bh;
    ///mass
    Double_t M_bh;
    ///mean accretion rate,metallicty
    Double_t acc_bh;
    //@}
#endif

    PropData(){
        num=gNFOF=0;
        gmass=gsize=gRmbp=gmaxvel=gRmaxvel=gRvir=gR200m=gR200c=gRhalfmass=Efrac=Pot=T=0.;
        gMFOF=0;
        gcm[0]=gcm[1]=gcm[2]=gcmvel[0]=gcmvel[1]=gcmvel[2]=0.;
        gJ[0]=gJ[1]=gJ[2]=0;
        gveldisp=Matrix(0.);
        gq=gs=1.0;
        Krot=0.;
#ifdef GASON
        n_gas=M_gas=Efrac_gas=0;
        cm_gas[0]=cm_gas[1]=cm_gas[2]=cmvel_gas[0]=cmvel_gas[1]=cmvel_gas[2]=0.;
        L_gas[0]=L_gas[1]=L_gas[2]=0;
        q_gas=s_gas=1.0;
        Rhalfmass_gas=0;
        eigvec_gas=Matrix(1,0,0,0,1,0,0,0,1);
        Temp_gas=Z_gas=SFR_gas=0.0;
        veldisp_gas=Matrix(0.);
        Krot_gas=T_gas=Pot_gas=0;
#endif
#ifdef STARON
        n_star=M_star=Efrac_star=0;
        cm_star[0]=cm_star[1]=cm_star[2]=cmvel_star[0]=cmvel_star[1]=cmvel_star[2]=0.;
        L_star[0]=L_star[1]=L_star[2]=0;
        q_star=s_star=1.0;
        Rhalfmass_star=0;
        eigvec_star=Matrix(1,0,0,0,1,0,0,0,1);
        t_star=Z_star=0.;
        veldisp_star=Matrix(0.);
        Krot_star=T_star=Pot_star=0;
#endif
#ifdef BHON
        n_bh=M_bh=0;
        acc_bh=0;
#endif
    }
    ///equals operator, useful if want inclusive information before substructure search
    PropData& operator=(const PropData &p){
        num=p.num;
        gcm=p.gcm;gcmvel=p.gcmvel;
        gpos=p.gpos;gvel=p.gvel;
        gmass=p.gmass;gsize=p.gsize;
        gMvir=p.gMvir;gRvir=p.gRvir;gRmbp=p.gRmbp;
        gmaxvel=gmaxvel=p.gmaxvel;gRmaxvel=p.gRmaxvel;gMmaxvel=p.gMmaxvel;
        gM200c=p.gM200c;gR200c=p.gR200c;
        gM200m=p.gM200m;gR200m=p.gR200m;
        return *this;
    }
    
    ///converts the properties data into comoving little h values
    ///so masses, positions have little h values and positions are comoving
    void ConverttoComove(Options &opt){
        gcm=gcm*opt.h/opt.a;
        gpos=gpos*opt.h/opt.a;
        gmass*=opt.h;
        gMvir*=opt.h;
        gM200c*=opt.h;
        gM200m*=opt.h;
        gMFOF*=opt.h;
        gsize*=opt.h/opt.a;
        gRmbp*=opt.h/opt.a;
        gRmaxvel*=opt.h/opt.a;
        gRvir*=opt.h/opt.a;
        gR200c*=opt.h/opt.a;
        gR200m*=opt.h/opt.a;
        gJ=gJ*opt.h*opt.h/opt.a;
#ifdef GASON
        M_gas*=opt.h;
        cm_gas=cm_gas*opt.h/opt.a;
        Rhalfmass_gas*=opt.h/opt.a;
        L_gas=L_gas*opt.h*opt.h/opt.a;
#endif
#ifdef STARON
        M_star*=opt.h;
        cm_star=cm_star*opt.h/opt.a;
        Rhalfmass_star*=opt.h/opt.a;
        L_star=L_star*opt.h*opt.h/opt.a;
#endif
#ifdef BHON
        M_bh*=opt.h;
#endif
    }

    void WriteBinary(fstream &Fout){
        long unsigned idval;
        unsigned int ival;
        double val, val3[3],val9[9];
        idval=haloid;
        Fout.write((char*)&idval,sizeof(idval));
        idval=ibound;
        Fout.write((char*)&idval,sizeof(idval));
        idval=hostid;
        Fout.write((char*)&idval,sizeof(idval));
        idval=numsubs;
        Fout.write((char*)&idval,sizeof(idval));
        idval=num;
        Fout.write((char*)&idval,sizeof(idval));
        
        val=gMvir;
        Fout.write((char*)&val,sizeof(val));
 
        for (int k=0;k<3;k++) val3[k]=gcm[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=gpos[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=gcmvel[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=gvel[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=gmass;
        Fout.write((char*)&val,sizeof(val));
        val=gMFOF;
        Fout.write((char*)&val,sizeof(val));
        val=gM200m;
        Fout.write((char*)&val,sizeof(val));
        val=gM200c;
        Fout.write((char*)&val,sizeof(val));
        val=gMvir;
        Fout.write((char*)&val,sizeof(val));

        val=Efrac;
        Fout.write((char*)&val,sizeof(val));

        val=gRvir;
        Fout.write((char*)&val,sizeof(val));
        val=gsize;
        Fout.write((char*)&val,sizeof(val));
        val=gR200m;
        Fout.write((char*)&val,sizeof(val));
        val=gR200c;
        Fout.write((char*)&val,sizeof(val));
        val=gRvir;
        Fout.write((char*)&val,sizeof(val));
        val=gRhalfmass;
        Fout.write((char*)&val,sizeof(val));
        val=gRmaxvel;
        Fout.write((char*)&val,sizeof(val));

        val=gmaxvel;
        Fout.write((char*)&val,sizeof(val));
        val=gsigma_v;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=gveldisp(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        val=glambda_B;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=gJ[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=gq;
        Fout.write((char*)&val,sizeof(val));
        val=gs;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=geigvec(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        val=cNFW;
        Fout.write((char*)&val,sizeof(val));
        val=Krot;
        Fout.write((char*)&val,sizeof(val));
        val=T;
        Fout.write((char*)&val,sizeof(val));
        val=Pot;
        Fout.write((char*)&val,sizeof(val));

#ifdef GASON
        idval=n_gas;
        Fout.write((char*)&idval,sizeof(idval));
        val=M_gas;
        Fout.write((char*)&val,sizeof(val));

        for (int k=0;k<3;k++) val3[k]=cm_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=cmvel_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=Efrac_gas;
        Fout.write((char*)&val,sizeof(val));

        val=Rhalfmass_gas;
        Fout.write((char*)&val,sizeof(val));

        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=veldisp_gas(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        for (int k=0;k<3;k++) val3[k]=L_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=q_gas;
        Fout.write((char*)&val,sizeof(val));
        val=s_gas;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=eigvec_gas(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        val=Krot_gas;
        Fout.write((char*)&val,sizeof(val));
        val=Temp_gas;
        Fout.write((char*)&val,sizeof(val));

#ifdef STARON
        val=Z_gas;
        Fout.write((char*)&val,sizeof(val));
        val=SFR_gas;
        Fout.write((char*)&val,sizeof(val));
#endif
#endif

#ifdef STARON
        idval=n_star;
        Fout.write((char*)&idval,sizeof(idval));
        val=M_star;
        Fout.write((char*)&val,sizeof(val));

        for (int k=0;k<3;k++) val3[k]=cm_star[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=cmvel_star[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=Efrac_star;
        Fout.write((char*)&val,sizeof(val));

        val=Rhalfmass_star;
        Fout.write((char*)&val,sizeof(val));

        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=veldisp_star(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        for (int k=0;k<3;k++) val3[k]=L_star[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=q_star;
        Fout.write((char*)&val,sizeof(val));
        val=s_star;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=eigvec_star(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        val=Krot_star;
        Fout.write((char*)&val,sizeof(val));
        val=t_star;
        Fout.write((char*)&val,sizeof(val));
        val=Z_star;
        Fout.write((char*)&val,sizeof(val));

#endif

        /*
        idbound=pdata[i].ibound;
        Fout.write((char*)&idbound,sizeof(long long));
        dvalue=pdata[i].gMvir;
        Fout.write((char*)&dvalue,sizeof(dvalue));
        ivalue=pdata[i].num;
        Fout.write((char*)&ivalue,sizeof(ivalue));
        for (int k=0;k<3;k++) ctemp[k]=pdata[i].gcm[k];
        Fout.write((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) ctemp[k]=pdata[i].gpos[k];
        Fout.write((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) ctemp[k]=pdata[i].gcmvel[k];
        Fout.write((char*)ctemp,sizeof(float)*3);
        for (int k=0;k<3;k++) ctemp[k]=pdata[i].gvel[k];
        Fout.write((char*)ctemp,sizeof(float)*3);
        value=pdata[i].gRvir;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].gRmbp;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].gRmaxvel;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].gmaxvel;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].gsigma_v;
        Fout.write((char*)&value,sizeof(value));
        for (int k=0;k<3;k++) ctemp[k]=pdata[i].gJ[k];
        Fout.write((char*)ctemp,sizeof(float)*3);
        value=pdata[i].gq;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].gs;
        Fout.write((char*)&value,sizeof(value));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) mtemp[k*3+n]=pdata[i].geigvec(k,n);
        Fout.write((char*)mtemp,sizeof(float)*9);
        //afterwards would be padding for 8 more floats (chars), add extra stuff like total mass, mass enclosed in Rmax, T, Pot, Efrac
        dvalue=pdata[i].gmass;
        Fout.write((char*)&dvalue,sizeof(dvalue));
        dvalue=pdata[i].gMmaxvel;
        Fout.write((char*)&dvalue,sizeof(dvalue));
        value=pdata[i].T;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].Pot;
        Fout.write((char*)&value,sizeof(value));
        value=pdata[i].Efrac;
        Fout.write((char*)&value,sizeof(value));
        value=0;
        for (int k=0;k<1;k++) Fout.write((char*)&value,sizeof(value));//for fact that 2*2*4+3*4
        */
    }
    
    void WriteAscii(fstream &Fout){
        Fout<<haloid<<" ";
        Fout<<ibound<<" ";
        Fout<<hostid<<" ";
        Fout<<numsubs<<" ";
        Fout<<num<<" ";
        Fout<<gMvir<<" ";
        for (int k=0;k<3;k++) Fout<<gcm[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gpos[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gcmvel[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gvel[k]<<" ";
        Fout<<gmass<<" ";
        Fout<<gMFOF<<" ";
        Fout<<gM200m<<" ";
        Fout<<gM200c<<" ";
        Fout<<gMvir<<" ";
        Fout<<Efrac<<" ";
        Fout<<gRvir<<" ";
        Fout<<gsize<<" ";
        Fout<<gR200m<<" ";
        Fout<<gR200c<<" ";
        Fout<<gRvir<<" ";
        Fout<<gRhalfmass<<" ";
        Fout<<gRmaxvel<<" ";
        Fout<<gmaxvel<<" ";
        Fout<<gsigma_v<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<gveldisp(k,n)<<" ";
        Fout<<glambda_B<<" ";
        for (int k=0;k<3;k++) Fout<<gJ[k]<<" ";
        Fout<<gq<<" ";
        Fout<<gs<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<geigvec(k,n)<<" ";
        Fout<<cNFW<<" ";
        Fout<<Krot<<" ";
        Fout<<T<<" ";
        Fout<<Pot<<" ";
#ifdef GASON
        Fout<<n_gas<<" ";
        Fout<<M_gas<<" ";
        for (int k=0;k<3;k++) Fout<<cm_gas[k]<<" ";
        for (int k=0;k<3;k++) Fout<<cmvel_gas[k]<<" ";
        Fout<<Efrac_gas<<" ";
        Fout<<Rhalfmass_gas<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<veldisp_gas(k,n)<<" ";
        for (int k=0;k<3;k++) Fout<<L_gas[k]<<" ";
        Fout<<q_gas<<" ";
        Fout<<s_gas<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<eigvec_gas(k,n)<<" ";
        Fout<<Krot_gas<<" ";
        Fout<<Temp_gas<<" ";
#ifdef STARON
        Fout<<Z_gas<<" ";
        Fout<<SFR_gas<<" ";
#endif
#endif

#ifdef STARON
        Fout<<n_star<<" ";
        Fout<<M_star<<" ";
        for (int k=0;k<3;k++) Fout<<cm_star[k]<<" ";
        for (int k=0;k<3;k++) Fout<<cmvel_star[k]<<" ";
        Fout<<Efrac_star<<" ";
        Fout<<Rhalfmass_star<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<veldisp_star(k,n)<<" ";
        for (int k=0;k<3;k++) Fout<<L_star[k]<<" ";
        Fout<<q_star<<" ";
        Fout<<s_star<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<eigvec_star(k,n)<<" ";
        Fout<<Krot_star<<" ";
        Fout<<t_star<<" ";
        Fout<<Z_star<<" ";
#endif
        Fout<<endl;
    }
#ifdef USEHDF
    void WriteHDF(fstream &Fout){};
#endif 
};

/*! Structures stores header info of the data writen by the \ref PropData data structure, 
    specifically the \ref PropData::WriteBinary, \ref PropData::WriteAscii, \ref PropData::WriteHDF routines
    Must ensure that these routines are all altered together so that the io makes sense.
*/
struct PropDataHeader{
    //list the header info
    vector<string> headerdatainfo;
    PropDataHeader(){
        headerdatainfo.push_back("ID");
        headerdatainfo.push_back("ID_mbp");
        headerdatainfo.push_back("hostHaloID");
        headerdatainfo.push_back("numSubStruct");
        headerdatainfo.push_back("npart");
        headerdatainfo.push_back("Mvir");
        headerdatainfo.push_back("Xc");
        headerdatainfo.push_back("Yc");
        headerdatainfo.push_back("Zc");
        headerdatainfo.push_back("Xcmbp");
        headerdatainfo.push_back("Ycmbp");
        headerdatainfo.push_back("Zcmbp");
        headerdatainfo.push_back("VXc");
        headerdatainfo.push_back("VYc");
        headerdatainfo.push_back("VZc");
        headerdatainfo.push_back("VXcmbp");
        headerdatainfo.push_back("VYcmbp");
        headerdatainfo.push_back("VZcmbp");
        headerdatainfo.push_back("Mass_tot");
        headerdatainfo.push_back("Mass_FOF");
        headerdatainfo.push_back("Mass_200mean");
        headerdatainfo.push_back("Mass_200crit");
        headerdatainfo.push_back("Mass_BN97");
        headerdatainfo.push_back("Efrac");
        headerdatainfo.push_back("Rvir");
        headerdatainfo.push_back("R_size");
        headerdatainfo.push_back("R_200mean");
        headerdatainfo.push_back("R_200crit");
        headerdatainfo.push_back("R_BN97");
        headerdatainfo.push_back("R_HalfMass");
        headerdatainfo.push_back("Rmax");
        headerdatainfo.push_back("Vmax");
        headerdatainfo.push_back("sigV");
        headerdatainfo.push_back("veldisp_xx");
        headerdatainfo.push_back("veldisp_xy");
        headerdatainfo.push_back("veldisp_xz");
        headerdatainfo.push_back("veldisp_yx");
        headerdatainfo.push_back("veldisp_yy");
        headerdatainfo.push_back("veldisp_yz");
        headerdatainfo.push_back("veldisp_zx");
        headerdatainfo.push_back("veldisp_zy");
        headerdatainfo.push_back("veldisp_zz");
        headerdatainfo.push_back("lambda_B");
        headerdatainfo.push_back("Lx");
        headerdatainfo.push_back("Ly");
        headerdatainfo.push_back("Lz");
        headerdatainfo.push_back("q");
        headerdatainfo.push_back("s");
        headerdatainfo.push_back("eig_xx");
        headerdatainfo.push_back("eig_xy");
        headerdatainfo.push_back("eig_xz");
        headerdatainfo.push_back("eig_yx");
        headerdatainfo.push_back("eig_yy");
        headerdatainfo.push_back("eig_yz");
        headerdatainfo.push_back("eig_zx");
        headerdatainfo.push_back("eig_zy");
        headerdatainfo.push_back("eig_zz");
        headerdatainfo.push_back("cNFW");
        headerdatainfo.push_back("Krot");
        headerdatainfo.push_back("Ekin");
        headerdatainfo.push_back("Epot");
#ifdef GASON
        headerdatainfo.push_back("n_gas");
        headerdatainfo.push_back("M_gas");
        headerdatainfo.push_back("Xc_gas");
        headerdatainfo.push_back("Yc_gas");
        headerdatainfo.push_back("Zc_gas");
        headerdatainfo.push_back("VXc_gas");
        headerdatainfo.push_back("VYc_gas");
        headerdatainfo.push_back("VZc_gas");
        headerdatainfo.push_back("Efrac_gas");
        headerdatainfo.push_back("R_HalfMass_gas");
        headerdatainfo.push_back("veldisp_xx_gas");
        headerdatainfo.push_back("veldisp_xy_gas");
        headerdatainfo.push_back("veldisp_xz_gas");
        headerdatainfo.push_back("veldisp_yx_gas");
        headerdatainfo.push_back("veldisp_yy_gas");
        headerdatainfo.push_back("veldisp_yz_gas");
        headerdatainfo.push_back("veldisp_zx_gas");
        headerdatainfo.push_back("veldisp_zy_gas");
        headerdatainfo.push_back("veldisp_zz_gas");
        headerdatainfo.push_back("Lx_gas");
        headerdatainfo.push_back("Ly_gas");
        headerdatainfo.push_back("Lz_gas");
        headerdatainfo.push_back("q_gas");
        headerdatainfo.push_back("s_gas");
        headerdatainfo.push_back("eig_xx_gas");
        headerdatainfo.push_back("eig_xy_gas");
        headerdatainfo.push_back("eig_xz_gas");
        headerdatainfo.push_back("eig_yx_gas");
        headerdatainfo.push_back("eig_yy_gas");
        headerdatainfo.push_back("eig_yz_gas");
        headerdatainfo.push_back("eig_zx_gas");
        headerdatainfo.push_back("eig_zy_gas");
        headerdatainfo.push_back("eig_zz_gas");
        headerdatainfo.push_back("Krot_gas");
        headerdatainfo.push_back("T_gas");
#ifdef STARON
        headerdatainfo.push_back("Zmet_gas");
        headerdatainfo.push_back("SFR_gas");
#endif
#endif
#ifdef STARON
        headerdatainfo.push_back("n_star");
        headerdatainfo.push_back("M_star");
        headerdatainfo.push_back("Xc_star");
        headerdatainfo.push_back("Yc_star");
        headerdatainfo.push_back("Zc_star");
        headerdatainfo.push_back("VXc_star");
        headerdatainfo.push_back("VYc_star");
        headerdatainfo.push_back("VZc_star");
        headerdatainfo.push_back("Efrac_star");
        headerdatainfo.push_back("R_HalfMass_star");
        headerdatainfo.push_back("veldisp_xx_star");
        headerdatainfo.push_back("veldisp_xy_star");
        headerdatainfo.push_back("veldisp_xz_star");
        headerdatainfo.push_back("veldisp_yx_star");
        headerdatainfo.push_back("veldisp_yy_star");
        headerdatainfo.push_back("veldisp_yz_star");
        headerdatainfo.push_back("veldisp_zx_star");
        headerdatainfo.push_back("veldisp_zy_star");
        headerdatainfo.push_back("veldisp_zz_star");
        headerdatainfo.push_back("Lx_star");
        headerdatainfo.push_back("Ly_star");
        headerdatainfo.push_back("Lz_star");
        headerdatainfo.push_back("q_star");
        headerdatainfo.push_back("s_star");
        headerdatainfo.push_back("eig_xx_star");
        headerdatainfo.push_back("eig_xy_star");
        headerdatainfo.push_back("eig_xz_star");
        headerdatainfo.push_back("eig_yx_star");
        headerdatainfo.push_back("eig_yy_star");
        headerdatainfo.push_back("eig_yz_star");
        headerdatainfo.push_back("eig_zx_star");
        headerdatainfo.push_back("eig_zy_star");
        headerdatainfo.push_back("eig_zz_star");
        headerdatainfo.push_back("Krot_star");
        headerdatainfo.push_back("tage_star");
        headerdatainfo.push_back("Zmet_star");
#endif
    }
};

/*! Structure used to keep track of a structure's parent structure
    note that level could be sim->halo->subhalo->subsubhalo
    or even sim->wall/void/filament->halo->substructure->subsubstructure
    here sim is stype=0,gid=0, and the other structure types are to be defined.
    The data structure is meant to be traversed from
        - level 0 structures ("field" objects)
        - level 0 pointer to nextlevel
        - nextlevel containing nginlevel objects
*/
struct StrucLevelData
{
    ///structure type and number in current level of hierarchy 
    Int_t stype,nsinlevel;
    //points to the the head pfof address of the group and parent
    Particle **Phead;
    Int_t **gidhead;
    //parent pointers point to the address of the parents gidhead and Phead
    Particle **Pparenthead;
    Int_t **gidparenthead;
    //add uber parent pointer (that is pointer to field halo)
    Int_t **giduberparenthead;
    ///allowing for multiple structure types at a given level in the hierarchy
    Int_t *stypeinlevel;
    StrucLevelData *nextlevel;
    StrucLevelData(Int_t numgroups=-1){
        if (numgroups<=0) {
            Phead=NULL;
            Pparenthead=NULL;
            gidhead=NULL;
            gidparenthead=NULL;
            giduberparenthead=NULL;
            nextlevel=NULL;
            stypeinlevel=NULL;
            nsinlevel=0;
        }
        else Allocate(numgroups);
    }
    //just allocate memory
    void Allocate(Int_t numgroups){
        nsinlevel=numgroups;
        Phead=new Particle*[numgroups+1];
        gidhead=new Int_t*[numgroups+1];
        stypeinlevel=new Int_t[numgroups+1];
        gidparenthead=new Int_t*[numgroups+1];
        giduberparenthead=new Int_t*[numgroups+1];
        nextlevel=NULL;
    }
    //initialize
    void Initialize(){
        for (Int_t i=1;i<=nsinlevel;i++) {gidhead[i]=NULL;gidparenthead[i]=NULL;giduberparenthead[i]=NULL;}
    }
};
extern StrucLevelData *psldata;

///if using MPI API
#ifdef USEMPI
#include <mpi.h>
///Includes external global variables used for MPI version of code
#include "mpivar.h"
#endif

#endif
