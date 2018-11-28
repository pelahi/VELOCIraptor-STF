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
#include <algorithm>
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

///if using ADIOS API
#ifdef USEADIOS
#include "adios.h"
#endif

//#include "swiftinterface.h"
//
//using namespace Swift;
using namespace std;
using namespace Math;
using namespace NBody;

//-- Structures and external variables

/// \defgroup PARTTYPES Particle types
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

/// \defgroup SEARCHTYPES Specify particle type to be searched, all, dm only, separate
//@{
#define  PSTALL 1
#define  PSTDARK 2
#define  PSTSTAR 3
#define  PSTGAS 4
#define  PSTBH 5
#define  PSTNOBH 6
//@}

/// \defgroup STRUCTURETYPES Specific structure type, allow for other types beside HALO
//@{
/// \todo note that here I have set background group type to a halo structure type but that can be changed
#define HALOSTYPE 10
#define HALOCORESTYPE 5
#define WALLSTYPE 1
#define VOIDSTYPE 2
#define FILAMENTSTYPE 3
#define BGTYPE 10
#define GROUPNOPARENT -1
#define FOF3DTYPE 7
#define FOF3DGROUP -2
//@}

/// \defgroup FOFTYPES FOF search types
//@{
//subsets made
///call \ref FOFStreamwithprob
#define  FOFSTPROB 1
///6D FOF search but only with outliers
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
//solely phase-space tensor core growth substructure search
#define  FOF6DCORE 6

///phase-space FOF but no subset produced
#define  FOFSTNOSUBSET 2
///no subsets made, just 6d (with each 6dfof search using 3d fof velocity dispersion,)
#define  FOF6DADAPTIVE 3
///6d fof but only use single velocity dispersion from largest 3d fof object
#define  FOF6D 4
///3d search
#define  FOF3D 5
//@}

/// \defgroup INTERATIVESEARCHPARAMS for iterative subsubstructure search
//@{
/// this is minimum particle number size for a subsearch to proceed whereby substructure split up into CELLSPLITNUM new cells
#define  MINCELLSIZE 100
#define  CELLSPLITNUM 8
#define  MINSUBSIZE MINCELLSIZE*CELLSPLITNUM
#define  MAXSUBLEVEL 8
/// maximum fraction a cell can take of a halo
#define  MAXCELLFRACTION 0.1
//@}

///\defgroup GRIDTYPES Type of Grid structures
//@{
#define  PHYSENGRID 1
#define  PHASEENGRID 2
#define  PHYSGRID 3
//@}

/// \name Max number of neighbouring cells used in interpolation of background velocity field.
//@{

//if cells were cubes this would be all the neighbours that full enclose the cube, 6 faces+20 diagonals
//using daganoals may not be ideal. Furthermore, code using adaptive grid that if effectively produces
//cell that are rectangular prisms. Furthermore, not all cells will share a boundary with another cell.
//So just consider "faces".
#define MAXNGRID 6

//@}

///\defgroup INPUTTYPES defining types of input
//@{
#define  NUMINPUTS 5
#define  IOGADGET 1
#define  IOHDF 2
#define  IOTIPSY 3
#define  IORAMSES 4
#define  IONCHILADA 5
//@}


///\defgroup OUTPUTTYPES defining format types of output
//@{
#define OUTASCII 0
#define OUTBINARY 1
#define OUTHDF 2
#define OUTADIOS 3
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

/// \defgroup OMPLIMS For determining whether loop contains enough for openm to be worthwhile.
//@{
#define ompsearchnum 50000
#define ompunbindnum 1000
#define ompperiodnum 50000
#define omppropnum 50000
//@}

/// \defgroup PROPLIMS Particle limits for calculating properties
//@{
#define propmincmnum 10
//@}

///\name halo id modifers used with current snapshot value to make temporally unique halo identifiers
#ifdef LONGINT
#define HALOIDSNVAL 1000000000000L
#else
#define HALOIDSNVAL 1000000
#endif

///\defgroup GASPARAMS Useful constants for gas
//@{
///mass of helium relative to hydrogen
#define M_HetoM_H 4.0026
//@}


/// Structure stores unbinding information
struct UnbindInfo
{
    ///\name flag whether unbind groups, keep bg potential when unbinding, type of unbinding and reference frame
    //@{
    int unbindflag,bgpot,unbindtype,cmvelreftype;
    //@}
    ///boolean as to whether code calculate potentials or potentials are externally provided
    bool icalculatepotential;
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
        icalculatepotential=true;
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

/* Structure to hold the location of a top-level cell. */
struct cell_loc {

    /* Coordinates x,y,z */
    double loc[3];

};

/// Options structure stores useful variables that have user determined values which are altered by \ref GetArgs in \ref ui.cxx
struct Options
{
    ///\name filenames
    //@{
    char *fname,*outname,*smname,*pname,*gname;
    char *ramsessnapname;
    //@}
    ///input format
    int inputtype;
    ///number of snapshots
    int num_files,snum;
    ///if parallel reading, number of files read in parallel
    int nsnapread;
    ///for output, specify the formats, ie. many separate files
    int iseparatefiles;
    ///for output specify the format HDF, binary or ascii \ref OUTHDF, \ref OUTBINARY, \ref OUTASCII
    int ibinaryout;
    ///for extended output allowing extraction of particles
    int iextendedoutput;
    /// output extra fields in halo properties
    int iextrahalooutput;
    /// calculate and output extra gas fields
    int iextragasoutput;
    /// calculate and output extra star fields
    int iextrastaroutput;
    /// calculate and output extra star fields
    int iextrabhoutput;

    ///disable particle id related output like fof.grp or catalog_group data. Useful if just want halo properties
    ///and not interested in tracking. Code writes halo properties catalog and exits.
    int inoidoutput;
    ///return propery data in in comoving little h units instead of standard physical units
    int icomoveunit;
    /// input is a cosmological simulation so can use box sizes, cosmological parameters, etc to set scales
    int icosmologicalin;
    /// input buffer size when reading data
    long int inputbufsize;
    /// mpi paritcle buffer size when sending input particle information
    long int mpiparticletotbufsize,mpiparticlebufsize;
    /// mpi factor by which to multiple the memory allocated, ie: buffer region
    /// to reduce likelihood of having to expand/allocate new memory
    Double_t mpipartfac;

    ///\name length,m,v,grav conversion units
    //@{
    Double_t L, M, U, V, G;
    Double_t lengthtokpc, velocitytokms, masstosolarmass;
    //@}
    ///period (comove)
    Double_t p;
    ///\name scale factor, Hubunit, h, cosmology, virial density. These are used if linking lengths are scaled or trying to define virlevel using the cosmology
    //@{
    Double_t a,H,h, Omega_m, Omega_b, Omega_cdm, Omega_Lambda, w_de, rhobg, virlevel, virBN98;
    int comove;
    /// to store the internal code unit to kpc and the distance^2 of 30 kpc, and 50 kpc
    Double_t lengthtokpc30pow2, lengthtokpc50pow2;
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
    ///\name parameters related to 3DFOF search & subsequent 6DFOF search
    //@{
    Double_t ellhalo6dxfac;
    Double_t ellhalo6dvfac;
    int iKeepFOF;
    Int_t num3dfof;
    //@}
    //@{
    ///\name factors used to check for halo mergers, large background substructures and store the velocity scale when searching for associated baryon substructures
    //@{
    Double_t HaloMergerSize,HaloMergerRatio,HaloSigmaV,HaloVelDispScale,HaloLocalSigmaV;
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

    ///if during substructure search, want to also search for larger substructures
    //using the more time consuming local velocity calculations (partly in lieu of using the faster core search)
    int iLargerCellSearch;

    ///\name extra stuff for halo merger check and identification of multiple halo core and flag for fully adaptive linking length using number density of candidate objects
    //@{
    /// run halo core search for mergers
    int iHaloCoreSearch;
    ///maximum sublevel at which we search for phase-space cores
    int maxnlevelcoresearch;
    ///parameters associated with phase-space search for cores of mergers
    Double_t halocorexfac, halocorevfac, halocorenfac, halocoresigmafac;
    ///x and v space linking lengths calculated for each object
    int iAdaptiveCoreLinking;
    ///use phase-space tensor core assignment
    int iPhaseCoreGrowth;
    ///number of iterations
    int halocorenumloops;
    ///factor by which one multiples the configuration space dispersion when looping for cores
    Double_t halocorexfaciter;
    ///factor by which one multiples the velocity space dispersion when looping for cores
    Double_t halocorevfaciter;
    ///factor by which one multiples the min num when looping for cores
    Double_t halocorenumfaciter;
    ///factor by which a core must be seperated from main core in phase-space in sigma units
    Double_t halocorephasedistsig;
    //@}
    ///for storing a snapshot value to make halo ids unique across snapshots
    long long snapshotvalue;

    ///\name for reading gadget info with lots of extra sph, star and bh blocks
    //@{
    int gnsphblocks,gnstarblocks,gnbhblocks;
    //@}

    /// \name Extra HDF flags indicating the existence of extra baryonic/dm particle types
    //@{
    /// input naming convention
    int ihdfnameconvention;
    /// input contains star particles
    int iusestarparticles;
    /// input contains black hole/sink particles
    int iusesinkparticles;
    /// input contains wind particles
    int iusewindparticles;
    /// input contains tracer particles
    int iusetracerparticles;
    /// input contains extra dark type particles
    int iuseextradarkparticles;
    //@}

    /// if want full spherical overdensity, factor by which size is multiplied to get
    ///bucket of particles
    Double_t SphericalOverdensitySeachFac;
    ///if want to the particle IDs that are within the SO overdensity of a halo
    int iSphericalOverdensityPartList;
    /// \name Extra variables to store information useful in zoom simluations
    //@{
    /// store the lowest dark matter particle mass
    Double_t zoomlowmassdm;
    //@}

    ///\name extra runtime flags
    //@{
    ///scale lengths. Useful if searching single halo system and which to automatically scale linking lengths
    int iScaleLengths;

    // Swift simulation information
    //Swift::siminfo swiftsiminfo;

    double spacedimension[3];

    /* Number of top-level cells. */
    int numcells;

    /* Number of top-level cells in each dimension. */
    int numcellsperdim;

    /* Locations of top-level cells. */
    cell_loc *cellloc;

    /*! Top-level cell width. */
    double cellwidth[3];

    /*! Inverse of the top-level cell width. */
    double icellwidth[3];

    /*! Holds the node ID of each top-level cell. */
    const int *cellnodeids;

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
        ellhalo6dxfac=1.0;
        ellhalo6dvfac=1.25;

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
        iKeepFOF=0;

        iLargerCellSearch=0;

        iHaloCoreSearch=0;
        iAdaptiveCoreLinking=0;
        iPhaseCoreGrowth=1;
        maxnlevelcoresearch=5;
        halocorexfac=0.5;
        halocorevfac=2.0;
        halocorenfac=0.1;
        halocoresigmafac=2.0;
        halocorenumloops=3;
        halocorexfaciter=0.75;
        halocorevfaciter=0.75;
        halocorenumfaciter=1.0;
        halocorephasedistsig=2.0;

        iverbose=0;
        iwritefof=0;
        iseparatefiles=0;
        ibinaryout=0;
        iextrahalooutput=0;
        iextendedoutput=0;
        inoidoutput=0;
        icomoveunit=0;
        icosmologicalin=1;

        iextragasoutput=0;
        iextrastaroutput=0;

        iusestarparticles=1;
        iusesinkparticles=1;
        iusewindparticles=0;
        iusetracerparticles=0;
#ifdef HIGHRES
        iuseextradarkparticles=1;
#else
        iuseextradarkparticles=0;
#endif

        snapshotvalue=0;

        gnsphblocks=4;
        gnstarblocks=2;
        gnbhblocks=2;

        iScaleLengths=0;

        inputbufsize=100000;

        mpiparticletotbufsize=-1;
        mpiparticlebufsize=-1;

        lengthtokpc=-1.0;
        velocitytokms=-1.0;
        masstosolarmass=-1.0;

        lengthtokpc30pow2=30.0*30.0;
        lengthtokpc30pow2=50.0*50.0;

        SphericalOverdensitySeachFac=1.25;
        iSphericalOverdensityPartList=0;

        mpipartfac=0.1;
#if USEHDF
        ihdfnameconvention=0;
#endif
    }
};

struct ConfigInfo{
    //list the name of the info
    vector<string> nameinfo;
    //vector<float> datainfo;
    vector<string> datainfo;
    //vector<int> datatype;
    ConfigInfo(Options &opt){
        //if compiler is super old and does not have at least std 11 implementation to_string does not exist
#ifndef OLDCCOMPILER
        //general search operations
        nameinfo.push_back("Particle_search_type");
        datainfo.push_back(to_string(opt.partsearchtype));
        nameinfo.push_back("FoF_search_type");
        datainfo.push_back(to_string(opt.foftype));
        nameinfo.push_back("FoF_Field_search_type");
        datainfo.push_back(to_string(opt.fofbgtype));
        nameinfo.push_back("Search_for_substructure");
        datainfo.push_back(to_string(opt.iSubSearch));
        nameinfo.push_back("Keep_FOF");
        datainfo.push_back(to_string(opt.iKeepFOF));
        nameinfo.push_back("Iterative_searchflag");
        datainfo.push_back(to_string(opt.iiterflag));
        nameinfo.push_back("Unbind_flag");
        datainfo.push_back(to_string(opt.uinfo.unbindflag));
        nameinfo.push_back("Baryon_searchflag");
        datainfo.push_back(to_string(opt.iBaryonSearch));
        nameinfo.push_back("CMrefadjustsubsearch_flag");
        datainfo.push_back(to_string(opt.icmrefadjust));
        nameinfo.push_back("Halo_core_search");
        datainfo.push_back(to_string(opt.iHaloCoreSearch));
        nameinfo.push_back("Use_adaptive_core_search");
        datainfo.push_back(to_string(opt.iAdaptiveCoreLinking));
        nameinfo.push_back("Use_phase_tensor_core_growth");
        datainfo.push_back(to_string(opt.iPhaseCoreGrowth));

        //local field parameters
        nameinfo.push_back("Cell_fraction");
        datainfo.push_back(to_string(opt.Ncellfac));
        nameinfo.push_back("Grid_type");
        datainfo.push_back(to_string(opt.gridtype));
        nameinfo.push_back("Nsearch_velocity");
        datainfo.push_back(to_string(opt.Nvel));
        nameinfo.push_back("Nsearch_physical");
        datainfo.push_back(to_string(opt.Nsearch));

        //substructure search parameters
        nameinfo.push_back("Outlier_threshold");
        datainfo.push_back(to_string(opt.ellthreshold));
        nameinfo.push_back("Significance_level");
        datainfo.push_back(to_string(opt.siglevel));
        nameinfo.push_back("Velocity_ratio");
        datainfo.push_back(to_string(opt.Vratio));
        nameinfo.push_back("Velocity_opening_angle");
        datainfo.push_back(to_string(opt.thetaopen));
        nameinfo.push_back("Physical_linking_length");
        datainfo.push_back(to_string(opt.ellphys));
        nameinfo.push_back("Velocity_linking_length");
        datainfo.push_back(to_string(opt.ellvel));
        nameinfo.push_back("Minimum_size");
        datainfo.push_back(to_string(opt.MinSize));

        //field object specific searches
        nameinfo.push_back("Minimum_halo_size");
        datainfo.push_back(to_string(opt.HaloMinSize));
        nameinfo.push_back("Halo_linking_length_factor");
        datainfo.push_back(to_string(opt.ellhalophysfac));
        nameinfo.push_back("Halo_velocity_linking_length_factor");
        datainfo.push_back(to_string(opt.ellhalovelfac));

        //specific to 6DFOF field search
        nameinfo.push_back("Halo_6D_linking_length_factor");
        datainfo.push_back(to_string(opt.ellhalo6dxfac));
        nameinfo.push_back("Halo_6D_vel_linking_length_factor");
        datainfo.push_back(to_string(opt.ellhalo6dvfac));

        //specific search for 6d fof core searches
        nameinfo.push_back("Halo_core_ellx_fac");
        datainfo.push_back(to_string(opt.halocorexfac));
        nameinfo.push_back("Halo_core_ellv_fac");
        datainfo.push_back(to_string(opt.halocorevfac));
        nameinfo.push_back("Halo_core_ncellfac");
        datainfo.push_back(to_string(opt.halocorenfac));
        nameinfo.push_back("Halo_core_adaptive_sigma_fac");
        datainfo.push_back(to_string(opt.halocoresigmafac));
        nameinfo.push_back("Halo_core_num_loops");
        datainfo.push_back(to_string(opt.halocorenumloops));
        nameinfo.push_back("Halo_core_loop_ellx_fac");
        datainfo.push_back(to_string(opt.halocorexfaciter));
        nameinfo.push_back("Halo_core_loop_ellv_fac");
        datainfo.push_back(to_string(opt.halocorevfaciter));
        nameinfo.push_back("Halo_core_loop_elln_fac");
        datainfo.push_back(to_string(opt.halocorenumfaciter));
        nameinfo.push_back("Halo_core_phase_significance");
        datainfo.push_back(to_string(opt.halocorephasedistsig));


        //for changing factors used in iterative search
        nameinfo.push_back("Iterative_threshold_factor");
        datainfo.push_back(to_string(opt.ellfac));
        nameinfo.push_back("Iterative_linking_length_factor");
        datainfo.push_back(to_string(opt.ellxfac));
        nameinfo.push_back("Iterative_Vratio_factor");
        datainfo.push_back(to_string(opt.vfac));
        nameinfo.push_back("Iterative_ThetaOp_factor");
        datainfo.push_back(to_string(opt.thetafac));

        //for changing effective resolution when rescaling linking lengh
        nameinfo.push_back("Effective_resolution");
        datainfo.push_back(to_string(opt.Neff));

        //for changing effective resolution when rescaling linking lengh
        nameinfo.push_back("Singlehalo_search");
        datainfo.push_back(to_string(opt.iSingleHalo));

        //units, cosmology
        nameinfo.push_back("Length_unit");
        datainfo.push_back(to_string(opt.L));
        nameinfo.push_back("Velocity_unit");
        datainfo.push_back(to_string(opt.V));
        nameinfo.push_back("Mass_unit");
        datainfo.push_back(to_string(opt.M));
        nameinfo.push_back("Hubble_unit");
        datainfo.push_back(to_string(opt.H));
        nameinfo.push_back("Gravity");
        datainfo.push_back(to_string(opt.G));
        nameinfo.push_back("Mass_value");
        datainfo.push_back(to_string(opt.MassValue));
        nameinfo.push_back("Length_unit_to_kpc");
        datainfo.push_back(to_string(opt.lengthtokpc));
        nameinfo.push_back("Velocity_to_kms");
        datainfo.push_back(to_string(opt.velocitytokms));
        nameinfo.push_back("Mass_to_solarmass");
        datainfo.push_back(to_string(opt.masstosolarmass));

        nameinfo.push_back("Period");
        datainfo.push_back(to_string(opt.p));
        nameinfo.push_back("Scale_factor");
        datainfo.push_back(to_string(opt.a));
        nameinfo.push_back("h_val");
        datainfo.push_back(to_string(opt.h));
        nameinfo.push_back("Omega_m");
        datainfo.push_back(to_string(opt.Omega_m));
        nameinfo.push_back("Omega_Lambda");
        datainfo.push_back(to_string(opt.Omega_Lambda));
        nameinfo.push_back("Critical_density");
        datainfo.push_back(to_string(opt.rhobg));
        nameinfo.push_back("Virial_density");
        datainfo.push_back(to_string(opt.virlevel));
        nameinfo.push_back("Omega_cdm");
        datainfo.push_back(to_string(opt.Omega_cdm));
        nameinfo.push_back("Omega_b");
        datainfo.push_back(to_string(opt.Omega_b));
        nameinfo.push_back("w_of_DE");
        datainfo.push_back(to_string(opt.w_de));

        //unbinding
        nameinfo.push_back("Softening_length");
        datainfo.push_back(to_string(opt.uinfo.eps));
        nameinfo.push_back("Allowed_kinetic_potential_ratio");
        datainfo.push_back(to_string(opt.uinfo.Eratio));
        nameinfo.push_back("Min_bound_mass_frac");
        datainfo.push_back(to_string(opt.uinfo.minEfrac));
        nameinfo.push_back("Bound_halos");
        datainfo.push_back(to_string(opt.iBoundHalos));
        nameinfo.push_back("Keep_background_potential");
        datainfo.push_back(to_string(opt.uinfo.bgpot));
        nameinfo.push_back("Kinetic_reference_frame_type");
        datainfo.push_back(to_string(opt.uinfo.cmvelreftype));
        nameinfo.push_back("Min_npot_ref");
        datainfo.push_back(to_string(opt.uinfo.Npotref));
        nameinfo.push_back("Frac_pot_ref");
        datainfo.push_back(to_string(opt.uinfo.fracpotref));
        nameinfo.push_back("Unbinding_type");
        datainfo.push_back(to_string(opt.uinfo.unbindtype));

        //other options
        nameinfo.push_back("Verbose");
        datainfo.push_back(to_string(opt.iverbose));
        nameinfo.push_back("Write_group_array_file");
        datainfo.push_back(to_string(opt.iwritefof));
        nameinfo.push_back("Snapshot_value");
        datainfo.push_back(to_string(opt.snapshotvalue));
        nameinfo.push_back("Inclusive_halo_masses");
        datainfo.push_back(to_string(opt.iInclusiveHalo));

        //io related
        nameinfo.push_back("Cosmological_input");
        datainfo.push_back(to_string(opt.icosmologicalin));
        nameinfo.push_back("Input_chunk_size");
        datainfo.push_back(to_string(opt.inputbufsize));
        nameinfo.push_back("MPI_particle_total_buf_size");
        datainfo.push_back(to_string(opt.mpiparticletotbufsize));
        nameinfo.push_back("Separate_output_files");
        datainfo.push_back(to_string(opt.iseparatefiles));
        nameinfo.push_back("Binary_output");
        datainfo.push_back(to_string(opt.ibinaryout));
        nameinfo.push_back("Comoving_units");
        datainfo.push_back(to_string(opt.icomoveunit));
        nameinfo.push_back("Extended_output");
        datainfo.push_back(to_string(opt.iextendedoutput));

        //gadget io related to extra info for sph, stars, bhs,
        nameinfo.push_back("NSPH_extra_blocks");
        datainfo.push_back(to_string(opt.gnsphblocks));
        nameinfo.push_back("NStar_extra_blocks");
        datainfo.push_back(to_string(opt.gnstarblocks));
        nameinfo.push_back("NBH_extra_blocks");
        datainfo.push_back(to_string(opt.gnbhblocks));

        //mpi related configuration
        nameinfo.push_back("MPI_part_allocation_fac");
        datainfo.push_back(to_string(opt.mpiparticletotbufsize));
#endif
    }
};

struct SimInfo{
    //list the name of the info
    vector<string> nameinfo;
    vector<string> datainfo;

    SimInfo(Options &opt){
        //if compiler is super old and does not have at least std 11 implementation to_string does not exist
#ifndef OLDCCOMPILER
        nameinfo.push_back("Cosmological_Sim");
        datainfo.push_back(to_string(opt.icosmologicalin));
        if (opt.icosmologicalin) {
            nameinfo.push_back("ScaleFactor");
            datainfo.push_back(to_string(opt.a));
            nameinfo.push_back("h_val");
            datainfo.push_back(to_string(opt.h));
            nameinfo.push_back("Omega_m");
            datainfo.push_back(to_string(opt.Omega_m));
            nameinfo.push_back("Omega_Lambda");
            datainfo.push_back(to_string(opt.Omega_Lambda));
            nameinfo.push_back("Omega_cdm");
            datainfo.push_back(to_string(opt.Omega_cdm));
            nameinfo.push_back("Omega_b");
            datainfo.push_back(to_string(opt.Omega_b));
            nameinfo.push_back("w_of_DE");
            datainfo.push_back(to_string(opt.w_de));
            nameinfo.push_back("Period");
            datainfo.push_back(to_string(opt.p));
            nameinfo.push_back("Hubble_unit");
            datainfo.push_back(to_string(opt.H));
        }
        else{
            nameinfo.push_back("Time");
            datainfo.push_back(to_string(opt.a));
            nameinfo.push_back("Period");
            datainfo.push_back(to_string(opt.p));
        }

        //units
        nameinfo.push_back("Length_unit");
        datainfo.push_back(to_string(opt.L));
        nameinfo.push_back("Velocity_unit");
        datainfo.push_back(to_string(opt.V));
        nameinfo.push_back("Mass_unit");
        datainfo.push_back(to_string(opt.M));
        nameinfo.push_back("Gravity");
        datainfo.push_back(to_string(opt.G));
#ifdef NOMASS
        nameinfo.push_back("Mass_value");
        datainfo.push_back(to_string(opt.MassValue));
#endif

#endif
    }
};


struct UnitInfo{
    //list the name of the info
    vector<string> nameinfo;
    vector<string> datainfo;

    UnitInfo(Options &opt){
        //if compiler is super old and does not have at least std 11 implementation to_string does not exist
#ifndef OLDCCOMPILER
        nameinfo.push_back("Cosmological_Sim");
        datainfo.push_back(to_string(opt.icosmologicalin));
        nameinfo.push_back("Comoving_or_Physical");
        datainfo.push_back(to_string(opt.icomoveunit));
        //units
        nameinfo.push_back("Length_unit_to_kpc");
        datainfo.push_back(to_string(opt.lengthtokpc));
        nameinfo.push_back("Velocity_unit_to_kms");
        datainfo.push_back(to_string(opt.velocitytokms));
        nameinfo.push_back("Mass_unit_to_solarmass");
        datainfo.push_back(to_string(opt.masstosolarmass));
#endif
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
    long long haloid,hostid,directhostid, hostfofid;
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
    Double_t gM200c,gR200c,gM200m,gR200m,gMFOF,gM500c,gR500c,gMBN98,gRBN98;
    //to store exclusive masses of halo ignoring substructure
    Double_t gMvir_excl,gRvir_excl,gM200c_excl,gR200c_excl,gM200m_excl,gR200m_excl,gMBN98_excl,gRBN98_excl;
    //@}
    ///\name physical properties for shape/mass distribution
    //@{
    ///axis ratios
    Double_t gq,gs;
    ///eigenvector
    Matrix geigvec;
    //@}
    ///\name physical properties for velocity
    //@{
    ///velocity dispersion
    Double_t gsigma_v;
    ///dispersion tensor
    Matrix gveldisp;
    //@}
    ///physical properties for dynamical state
    Double_t Efrac,Pot,T;
    ///physical properties for angular momentum
    Coordinate gJ;
    Coordinate gJ200m, gJ200c, gJBN98;
    ///physical properties for angular momentum exclusive
    Coordinate gJ200m_excl, gJ200c_excl, gJBN98_excl;
    ///Keep track of position of least unbound particle and most bound particle pid
    Int_t iunbound,ibound;
    ///Type of structure
    int stype;
    ///concentration (and related quantity used to calculate a concentration)
    Double_t cNFW, VmaxVvir2;
    ///Bullock & Peebles spin parameters
    Double_t glambda_B,glambda_P;
    ///measure of rotational support
    Double_t Krot;
    //@}

    ///\name halo properties within RVmax
    //@{
    Double_t RV_q,RV_s;
    Matrix RV_eigvec;
    Double_t RV_sigma_v;
    Matrix RV_veldisp;
    Coordinate RV_J;
    Double_t RV_lambda_B,RV_lambda_P;
    Double_t RV_Krot;
    //@}

#ifdef GASON
    ///\name gas specific quantities
    //@{
    ///number of particles
    int n_gas;
    ///mass
    Double_t M_gas, M_gas_rvmax,M_gas_30kpc,M_gas_50kpc, M_gas_500c;
    ///mass in spherical overdensities
    Double_t M_200crit_gas, M_200mean_gas, M_BN98_gas;
    ///mass in spherical overdensities inclusive of all masses
    Double_t M_200crit_excl_gas, M_200mean_excl_gas, M_BN98_excl_gas;
    ///pos/vel info
    Coordinate cm_gas,cmvel_gas;
    ///velocity/angular momentum info
    Double_t Krot_gas;
    Coordinate L_gas;
    ///physical properties for angular momentum (can be inclusive or exclusive )
    Coordinate L_200crit_gas, L_200mean_gas, L_BN98_gas;
    ///physical properties for angular momentum exclusiveto object
    Coordinate L_200crit_excl_gas, L_200mean_excl_gas, L_BN98_excl_gas;
    //dispersion
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
    Double_t M_star, M_star_rvmax, M_star_30kpc, M_star_50kpc, M_star_500c;
    ///mass in spherical overdensities
    Double_t M_200crit_star, M_200mean_star, M_BN98_star;
    ///mass in spherical overdensities inclusive of all masses
    Double_t M_200crit_excl_star, M_200mean_excl_star, M_BN98_excl_star;
    ///pos/vel info
    Coordinate cm_star,cmvel_star;
    ///velocity/angular momentum info
    Double_t Krot_star;
    Coordinate L_star;
    ///physical properties for angular momentum (can be inclusive or exclusive )
    Coordinate L_200crit_star, L_200mean_star, L_BN98_star;
    ///physical properties for angular momentum exclusiveto object
    Coordinate L_200crit_excl_star, L_200mean_excl_star, L_BN98_excl_star;
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

#ifdef HIGHRES
    ///\name low resolution interloper particle specific quantities
    //@{
    ///number of interloper low res particles
    int n_interloper;
    ///mass
    Double_t M_interloper;
    //@}
#endif

    PropData(){
        num=gNFOF=0;
        gmass=gsize=gRmbp=gmaxvel=gRmaxvel=gRvir=gR200m=gR200c=gRhalfmass=Efrac=Pot=T=0.;
        gMFOF=0;
        gM500c=gR500c=0;
        gMBN98=gRBN98=0;
        gcm[0]=gcm[1]=gcm[2]=gcmvel[0]=gcmvel[1]=gcmvel[2]=0.;
        gJ[0]=gJ[1]=gJ[2]=0;
        gJ200m[0]=gJ200m[1]=gJ200m[2]=0;
        gJ200c[0]=gJ200c[1]=gJ200c[2]=0;
        gJBN98[0]=gJBN98[1]=gJBN98[2]=0;
        gveldisp=Matrix(0.);
        gq=gs=1.0;
        Krot=0.;

        gM200m_excl=gM200c_excl=gMBN98_excl=0;
        gR200m_excl=gR200c_excl=gRBN98_excl=0;
        gJ200m_excl[0]=gJ200m_excl[1]=gJ200m_excl[2]=0;
        gJ200c_excl[0]=gJ200c_excl[1]=gJ200c_excl[2]=0;
        gJBN98_excl[0]=gJBN98_excl[1]=gJBN98_excl[2]=0;

        RV_sigma_v=0;
        RV_q=RV_s=1.;
        RV_J[0]=RV_J[1]=RV_J[2]=0;
        RV_veldisp=Matrix(0.);
        RV_eigvec=Matrix(0.);
        RV_lambda_B=RV_lambda_P=RV_Krot=0;

#ifdef GASON
        M_gas_rvmax=M_gas_30kpc=M_gas_50kpc=0;
        n_gas=M_gas=Efrac_gas=0;
        cm_gas[0]=cm_gas[1]=cm_gas[2]=cmvel_gas[0]=cmvel_gas[1]=cmvel_gas[2]=0.;
        L_gas[0]=L_gas[1]=L_gas[2]=0;
        q_gas=s_gas=1.0;
        Rhalfmass_gas=0;
        eigvec_gas=Matrix(1,0,0,0,1,0,0,0,1);
        Temp_gas=Z_gas=SFR_gas=0.0;
        veldisp_gas=Matrix(0.);
        Krot_gas=T_gas=Pot_gas=0;

        M_200mean_gas=M_200crit_gas=M_BN98_gas=0;
        M_200mean_excl_gas=M_200crit_excl_gas=M_BN98_excl_gas=0;
        L_200crit_gas[0]=L_200crit_gas[1]=L_200crit_gas[2]=0;
        L_200mean_gas[0]=L_200mean_gas[1]=L_200mean_gas[2]=0;
        L_BN98_gas[0]=L_BN98_gas[1]=L_BN98_gas[2]=0;
        L_200crit_excl_gas[0]=L_200crit_excl_gas[1]=L_200crit_excl_gas[2]=0;
        L_200mean_excl_gas[0]=L_200mean_excl_gas[1]=L_200mean_excl_gas[2]=0;
        L_BN98_excl_gas[0]=L_BN98_excl_gas[1]=L_BN98_excl_gas[2]=0;
#endif
#ifdef STARON
        M_star_rvmax=M_star_30kpc=M_star_50kpc=0;
        n_star=M_star=Efrac_star=0;
        cm_star[0]=cm_star[1]=cm_star[2]=cmvel_star[0]=cmvel_star[1]=cmvel_star[2]=0.;
        L_star[0]=L_star[1]=L_star[2]=0;
        q_star=s_star=1.0;
        Rhalfmass_star=0;
        eigvec_star=Matrix(1,0,0,0,1,0,0,0,1);
        t_star=Z_star=0.;
        veldisp_star=Matrix(0.);
        Krot_star=T_star=Pot_star=0;

        M_200mean_star=M_200crit_star=M_BN98_star=0;
        M_200mean_excl_star=M_200crit_excl_star=M_BN98_excl_star=0;
        L_200crit_star[0]=L_200crit_star[1]=L_200crit_star[2]=0;
        L_200mean_star[0]=L_200mean_star[1]=L_200mean_star[2]=0;
        L_BN98_star[0]=L_BN98_star[1]=L_BN98_star[2]=0;
        L_200crit_excl_star[0]=L_200crit_excl_star[1]=L_200crit_excl_star[2]=0;
        L_200mean_excl_star[0]=L_200mean_excl_star[1]=L_200mean_excl_star[2]=0;
        L_BN98_excl_star[0]=L_BN98_excl_star[1]=L_BN98_excl_star[2]=0;
#endif
#ifdef BHON
        n_bh=M_bh=0;
        acc_bh=0;
#endif
#ifdef HIGHRES
        n_interloper=M_interloper=0;
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
        gM500c=p.gM500c;gR500c=p.gR500c;
        gMBN98=p.gMBN98;gRBN98=p.gRBN98;
        gNFOF=p.gNFOF;
        gMFOF=p.gMFOF;

        gM200c_excl=p.gM200c_excl;gR200c_excl=p.gR200c_excl;
        gM200m_excl=p.gM200m_excl;gR200m_excl=p.gR200m_excl;
        gMBN98_excl=p.gMBN98_excl;gRBN98_excl=p.gRBN98_excl;
        gJ=p.gJ;
        gJ200c=p.gJ200c;
        gJ200m=p.gJ200m;
        gJBN98=p.gJBN98;
        gJ200c_excl=p.gJ200c_excl;
        gJ200m_excl=p.gJ200m_excl;
        gJBN98_excl=p.gJBN98_excl;

        ///expand to copy all the gas, star, bh, stuff
#ifdef GASON
        M_200mean_gas=p.M_200mean_gas;
        M_200crit_gas=p.M_200crit_gas;
        M_BN98_gas=p.M_BN98_gas;
        M_200mean_excl_gas=p.M_200mean_excl_gas;
        M_200crit_excl_gas=p.M_200crit_excl_gas;
        M_BN98_excl_gas=p.M_BN98_excl_gas;
        L_200mean_gas=p.L_200mean_gas;
        L_200crit_gas=p.L_200crit_gas;
        L_BN98_gas=p.L_BN98_gas;
        L_200mean_excl_gas=p.L_200mean_excl_gas;
        L_200crit_excl_gas=p.L_200crit_excl_gas;
        L_BN98_excl_gas=p.L_BN98_excl_gas;
#endif
#ifdef STARON
        M_200mean_star=p.M_200mean_star;
        M_200crit_star=p.M_200crit_star;
        M_BN98_star=p.M_BN98_star;
        M_200mean_excl_star=p.M_200mean_excl_star;
        M_200crit_excl_star=p.M_200crit_excl_star;
        M_BN98_excl_star=p.M_BN98_excl_star;
        L_200mean_star=p.L_200mean_star;
        L_200crit_star=p.L_200crit_star;
        L_BN98_star=p.L_BN98_star;
        L_200mean_excl_star=p.L_200mean_excl_star;
        L_200crit_excl_star=p.L_200crit_excl_star;
        L_BN98_excl_star=p.L_BN98_excl_star;
#endif
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
        gM500c*=opt.h;
        gMBN98*=opt.h;
        gMFOF*=opt.h;
        gsize*=opt.h/opt.a;
        gRmbp*=opt.h/opt.a;
        gRmaxvel*=opt.h/opt.a;
        gRvir*=opt.h/opt.a;
        gR200c*=opt.h/opt.a;
        gR200m*=opt.h/opt.a;
        gR500c*=opt.h/opt.a;
        gRBN98*=opt.h/opt.a;
        gJ=gJ*opt.h*opt.h/opt.a;
        gJ200m=gJ200m*opt.h*opt.h/opt.a;
        gJ200c=gJ200c*opt.h*opt.h/opt.a;
        gJBN98=gJBN98*opt.h*opt.h/opt.a;
        RV_J=RV_J*opt.h*opt.h/opt.a;

        if (opt.iextrahalooutput) {
            gM200c_excl*=opt.h;
            gM200m_excl*=opt.h;
            gMBN98_excl*=opt.h;
            gR200c_excl*=opt.h/opt.a;
            gR200m_excl*=opt.h/opt.a;
            gRBN98_excl*=opt.h/opt.a;
            gJ200m_excl=gJ200m_excl*opt.h*opt.h/opt.a;
            gJ200c_excl=gJ200c_excl*opt.h*opt.h/opt.a;
            gJBN98_excl=gJBN98_excl*opt.h*opt.h/opt.a;
        }

#ifdef GASON
        M_gas*=opt.h;
        M_gas_rvmax*=opt.h;
        M_gas_30kpc*=opt.h;
        M_gas_50kpc*=opt.h;
        M_gas_500c*=opt.h;

        cm_gas=cm_gas*opt.h/opt.a;
        Rhalfmass_gas*=opt.h/opt.a;
        L_gas=L_gas*opt.h*opt.h/opt.a;

        if (opt.iextragasoutput) {
            M_200mean_gas*=opt.h;
            M_200crit_gas*=opt.h;
            M_BN98_gas*=opt.h;
            M_200mean_excl_gas*=opt.h;
            M_200crit_excl_gas*=opt.h;
            M_BN98_excl_gas*=opt.h;
            L_200crit_gas=L_200crit_gas*opt.h*opt.h/opt.a;
            L_200mean_gas=L_200mean_gas*opt.h*opt.h/opt.a;
            L_BN98_gas=L_BN98_gas*opt.h*opt.h/opt.a;
            L_200crit_excl_gas=L_200crit_excl_gas*opt.h*opt.h/opt.a;
            L_200mean_excl_gas=L_200mean_excl_gas*opt.h*opt.h/opt.a;
            L_BN98_excl_gas=L_BN98_excl_gas*opt.h*opt.h/opt.a;
        }

#endif
#ifdef STARON
        M_star*=opt.h;
        M_star_rvmax*=opt.h;
        M_star_30kpc*=opt.h;
        M_star_50kpc*=opt.h;
        M_star_500c*=opt.h;
        cm_star=cm_star*opt.h/opt.a;
        Rhalfmass_star*=opt.h/opt.a;
        L_star=L_star*opt.h*opt.h/opt.a;
#endif
#ifdef BHON
        M_bh*=opt.h;
#endif
#ifdef HIGHRES
        M_interloper*=opt.h;
#endif
    }

    ///write (append) the properties data to an already open binary file
    void WriteBinary(fstream &Fout, Options&opt){
        long long lval;
        long unsigned idval;
        unsigned int ival;
        double val, val3[3],val9[9];
        idval=haloid;
        Fout.write((char*)&idval,sizeof(idval));
        lval=ibound;
        Fout.write((char*)&lval,sizeof(idval));
        lval=hostid;
        Fout.write((char*)&lval,sizeof(idval));
        idval=numsubs;
        Fout.write((char*)&idval,sizeof(idval));
        idval=num;
        Fout.write((char*)&idval,sizeof(idval));
        ival=stype;
        Fout.write((char*)&ival,sizeof(ival));
        if (opt.iKeepFOF==1) {
            idval=directhostid;
            Fout.write((char*)&idval,sizeof(idval));
            idval=hostfofid;
            Fout.write((char*)&idval,sizeof(idval));
        }

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
        val=gMBN98;
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
        val=gRBN98;
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

        val=RV_sigma_v;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=RV_veldisp(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        val=RV_lambda_B;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=RV_J[k];
        Fout.write((char*)val3,sizeof(val)*3);

        val=RV_q;
        Fout.write((char*)&val,sizeof(val));
        val=RV_s;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=RV_eigvec(k,n);
        Fout.write((char*)val9,sizeof(val)*9);

        if (opt.iextrahalooutput) {
            for (int k=0;k<3;k++) val3[k]=gJ200m[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=gJ200c[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=gJBN98[k];
            Fout.write((char*)val3,sizeof(val)*3);
            if (opt.iInclusiveHalo>0) {
                val=gM200m_excl;
                Fout.write((char*)&val,sizeof(val));
                val=gM200c_excl;
                Fout.write((char*)&val,sizeof(val));
                val=gMBN98_excl;
                Fout.write((char*)&val,sizeof(val));

                val=gR200m_excl;
                Fout.write((char*)&val,sizeof(val));
                val=gR200c_excl;
                Fout.write((char*)&val,sizeof(val));
                val=gRBN98_excl;
                Fout.write((char*)&val,sizeof(val));

                for (int k=0;k<3;k++) val3[k]=gJ200m_excl[k];
                Fout.write((char*)val3,sizeof(val)*3);
                for (int k=0;k<3;k++) val3[k]=gJ200c_excl[k];
                Fout.write((char*)val3,sizeof(val)*3);
                for (int k=0;k<3;k++) val3[k]=gJBN98_excl[k];
                Fout.write((char*)val3,sizeof(val)*3);
            }
        }
#ifdef GASON
        idval=n_gas;
        Fout.write((char*)&idval,sizeof(idval));
        val=M_gas;
        Fout.write((char*)&val,sizeof(val));
        val=M_gas_rvmax;
        Fout.write((char*)&val,sizeof(val));
        val=M_gas_30kpc;
        Fout.write((char*)&val,sizeof(val));
        val=M_gas_500c;
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

    if (opt.iextragasoutput) {
        val=M_200mean_gas;
        Fout.write((char*)&val,sizeof(val));
        val=M_200crit_gas;
        Fout.write((char*)&val,sizeof(val));
        val=M_BN98_gas;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=L_200mean_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_200crit_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_BN98_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
        if (opt.iInclusiveHalo>0) {
            val=M_200mean_excl_gas;
            Fout.write((char*)&val,sizeof(val));
            val=M_200crit_excl_gas;
            Fout.write((char*)&val,sizeof(val));
            val=M_BN98_excl_gas;
            Fout.write((char*)&val,sizeof(val));
            for (int k=0;k<3;k++) val3[k]=L_200mean_excl_gas[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_200crit_excl_gas[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_BN98_excl_gas[k];
            Fout.write((char*)val3,sizeof(val)*3);
        }
    }
#endif

#ifdef STARON
        idval=n_star;
        Fout.write((char*)&idval,sizeof(idval));
        val=M_star;
        Fout.write((char*)&val,sizeof(val));
        val=M_star_rvmax;
        Fout.write((char*)&val,sizeof(val));
        val=M_star_30kpc;
        Fout.write((char*)&val,sizeof(val));
        val=M_star_500c;
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

        if (opt.iextragasoutput) {
            val=M_200mean_star;
            Fout.write((char*)&val,sizeof(val));
            val=M_200crit_star;
            Fout.write((char*)&val,sizeof(val));
            val=M_BN98_star;
            Fout.write((char*)&val,sizeof(val));
            for (int k=0;k<3;k++) val3[k]=L_200mean_star[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_200crit_star[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_BN98_star[k];
            Fout.write((char*)val3,sizeof(val)*3);
            if (opt.iInclusiveHalo>0) {
                val=M_200mean_excl_star;
                Fout.write((char*)&val,sizeof(val));
                val=M_200crit_excl_star;
                Fout.write((char*)&val,sizeof(val));
                val=M_BN98_excl_star;
                Fout.write((char*)&val,sizeof(val));
                for (int k=0;k<3;k++) val3[k]=L_200mean_excl_star[k];
                Fout.write((char*)val3,sizeof(val)*3);
                for (int k=0;k<3;k++) val3[k]=L_200crit_excl_star[k];
                Fout.write((char*)val3,sizeof(val)*3);
                for (int k=0;k<3;k++) val3[k]=L_BN98_excl_star[k];
                Fout.write((char*)val3,sizeof(val)*3);
            }
        }
#endif

#ifdef BHON
        idval=n_bh;
        Fout.write((char*)&idval,sizeof(idval));
        val=M_bh;
        Fout.write((char*)&val,sizeof(val));
#endif
#ifdef HIGHRES
        idval=n_interloper;
        Fout.write((char*)&idval,sizeof(idval));
        val=M_interloper;
        Fout.write((char*)&val,sizeof(val));
#endif
    }

    ///write (append) the properties data to an already open ascii file
    void WriteAscii(fstream &Fout, Options&opt){
        Fout<<haloid<<" ";
        Fout<<ibound<<" ";
        Fout<<hostid<<" ";
        Fout<<numsubs<<" ";
        Fout<<num<<" ";
        Fout<<stype<<" ";
        if (opt.iKeepFOF==1) {
            Fout<<directhostid<<" ";
            Fout<<hostfofid<<" ";
        }
        Fout<<gMvir<<" ";
        for (int k=0;k<3;k++) Fout<<gcm[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gpos[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gcmvel[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gvel[k]<<" ";
        Fout<<gmass<<" ";
        Fout<<gMFOF<<" ";
        Fout<<gM200m<<" ";
        Fout<<gM200c<<" ";
        Fout<<gMBN98<<" ";
        Fout<<Efrac<<" ";
        Fout<<gRvir<<" ";
        Fout<<gsize<<" ";
        Fout<<gR200m<<" ";
        Fout<<gR200c<<" ";
        Fout<<gRBN98<<" ";
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

        Fout<<RV_sigma_v<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<RV_veldisp(k,n)<<" ";
        Fout<<RV_lambda_B<<" ";
        for (int k=0;k<3;k++) Fout<<RV_J[k]<<" ";
        Fout<<RV_q<<" ";
        Fout<<RV_s<<" ";
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<RV_eigvec(k,n)<<" ";

        if (opt.iextrahalooutput) {
            for (int k=0;k<3;k++) Fout<<gJ200m[k]<<" ";
            for (int k=0;k<3;k++) Fout<<gJ200c[k]<<" ";
            for (int k=0;k<3;k++) Fout<<gJBN98[k]<<" ";
            if (opt.iInclusiveHalo>0) {
                Fout<<gM200m_excl<<" ";
                Fout<<gM200c_excl<<" ";
                Fout<<gMBN98_excl<<" ";
                Fout<<gR200m_excl<<" ";
                Fout<<gR200c_excl<<" ";
                Fout<<gRBN98_excl<<" ";
                for (int k=0;k<3;k++) Fout<<gJ200m_excl[k]<<" ";
                for (int k=0;k<3;k++) Fout<<gJ200c_excl[k]<<" ";
                for (int k=0;k<3;k++) Fout<<gJBN98_excl[k]<<" ";
            }
        }

#ifdef GASON
        Fout<<n_gas<<" ";
        Fout<<M_gas<<" ";
        Fout<<M_gas_rvmax<<" ";
        Fout<<M_gas_30kpc<<" ";
        //Fout<<M_gas_50kpc<<" ";
        Fout<<M_gas_500c<<" ";
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
    if (opt.iextragasoutput) {
        Fout<<M_200mean_gas<<" ";
        Fout<<M_200crit_gas<<" ";
        Fout<<M_BN98_gas<<" ";
        for (int k=0;k<3;k++) Fout<<L_200mean_gas[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_200crit_gas[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_BN98_gas[k]<<" ";
        if (opt.iInclusiveHalo>0) {
            Fout<<M_200mean_excl_gas<<" ";
            Fout<<M_200crit_excl_gas<<" ";
            Fout<<M_BN98_excl_gas<<" ";
            for (int k=0;k<3;k++) Fout<<L_200mean_excl_gas[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_200crit_excl_gas[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_BN98_excl_gas[k]<<" ";
        }
    }
#endif

#ifdef STARON
        Fout<<n_star<<" ";
        Fout<<M_star<<" ";
        Fout<<M_star_rvmax<<" ";
        Fout<<M_star_30kpc<<" ";
        //Fout<<M_star_50kpc<<" ";
        Fout<<M_star_500c<<" ";
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
        if (opt.iextragasoutput) {
            Fout<<M_200mean_star<<" ";
            Fout<<M_200crit_star<<" ";
            Fout<<M_BN98_star<<" ";
            for (int k=0;k<3;k++) Fout<<L_200mean_star[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_200crit_star[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_BN98_star[k]<<" ";
            if (opt.iInclusiveHalo>0) {
                Fout<<M_200mean_excl_star<<" ";
                Fout<<M_200crit_excl_star<<" ";
                Fout<<M_BN98_excl_star<<" ";
                for (int k=0;k<3;k++) Fout<<L_200mean_excl_star[k]<<" ";
                for (int k=0;k<3;k++) Fout<<L_200crit_excl_star[k]<<" ";
                for (int k=0;k<3;k++) Fout<<L_BN98_excl_star[k]<<" ";
            }
        }
#endif

#ifdef BHON
        Fout<<n_bh<<" ";
        Fout<<M_bh<<" ";
#endif
#ifdef HIGHRES
        Fout<<n_interloper<<" ";
        Fout<<M_interloper<<" ";
#endif
        Fout<<endl;
    }
#ifdef USEHDF
    ///write (append) the properties data to an already open hdf file
    void WriteHDF(H5File &Fhdf, DataSpace *&dataspaces, DataSet *&datasets, Options&opt){
    };
#endif
};

/*! Structures stores header info of the data writen by the \ref PropData data structure,
    specifically the \ref PropData::WriteBinary, \ref PropData::WriteAscii, \ref PropData::WriteHDF routines
    Must ensure that these routines are all altered together so that the io makes sense.
*/
struct PropDataHeader{
    //list the header info
    vector<string> headerdatainfo;
#ifdef USEHDF
    vector<PredType> predtypeinfo;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospredtypeinfo;
#endif
    PropDataHeader(Options&opt){
        int sizeval;
#ifdef USEHDF
        vector<PredType> desiredproprealtype;
        if (sizeof(Double_t)==sizeof(double)) desiredproprealtype.push_back(PredType::NATIVE_DOUBLE);
        else desiredproprealtype.push_back(PredType::NATIVE_FLOAT);
#endif
#ifdef USEADIOS
        vector<ADIOS_DATATYPES> desiredadiosproprealtype;
        if (sizeof(Double_t)==sizeof(double)) desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_double);
        else desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_real);
#endif

        headerdatainfo.push_back("ID");
        headerdatainfo.push_back("ID_mbp");
        headerdatainfo.push_back("hostHaloID");
        headerdatainfo.push_back("numSubStruct");
        headerdatainfo.push_back("npart");
        headerdatainfo.push_back("Structuretype");
        if (opt.iKeepFOF==1){
            headerdatainfo.push_back("hostDirectHaloID");
            headerdatainfo.push_back("hostFOFID");
        }

        //if using hdf, store the type
#ifdef USEHDF
        predtypeinfo.push_back(PredType::STD_U64LE);
        predtypeinfo.push_back(PredType::STD_I64LE);
        predtypeinfo.push_back(PredType::STD_I64LE);
        predtypeinfo.push_back(PredType::STD_U64LE);
        predtypeinfo.push_back(PredType::STD_U64LE);
        predtypeinfo.push_back(PredType::STD_I32LE);
        if (opt.iKeepFOF==1){
            predtypeinfo.push_back(PredType::STD_I64LE);
            predtypeinfo.push_back(PredType::STD_I64LE);
        }
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_integer);
        if (opt.iKeepFOF==1){
            adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
            adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
        }
#endif

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
        headerdatainfo.push_back("Mass_BN98");
        headerdatainfo.push_back("Efrac");
        headerdatainfo.push_back("Rvir");
        headerdatainfo.push_back("R_size");
        headerdatainfo.push_back("R_200mean");
        headerdatainfo.push_back("R_200crit");
        headerdatainfo.push_back("R_BN98");
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

        //some properties within RVmax
        headerdatainfo.push_back("RVmax_sigV");
        headerdatainfo.push_back("RVmax_veldisp_xx");
        headerdatainfo.push_back("RVmax_veldisp_xy");
        headerdatainfo.push_back("RVmax_veldisp_xz");
        headerdatainfo.push_back("RVmax_veldisp_yx");
        headerdatainfo.push_back("RVmax_veldisp_yy");
        headerdatainfo.push_back("RVmax_veldisp_yz");
        headerdatainfo.push_back("RVmax_veldisp_zx");
        headerdatainfo.push_back("RVmax_veldisp_zy");
        headerdatainfo.push_back("RVmax_veldisp_zz");
        headerdatainfo.push_back("RVmax_lambda_B");
        headerdatainfo.push_back("RVmax_Lx");
        headerdatainfo.push_back("RVmax_Ly");
        headerdatainfo.push_back("RVmax_Lz");
        headerdatainfo.push_back("RVmax_q");
        headerdatainfo.push_back("RVmax_s");
        headerdatainfo.push_back("RVmax_eig_xx");
        headerdatainfo.push_back("RVmax_eig_xy");
        headerdatainfo.push_back("RVmax_eig_xz");
        headerdatainfo.push_back("RVmax_eig_yx");
        headerdatainfo.push_back("RVmax_eig_yy");
        headerdatainfo.push_back("RVmax_eig_yz");
        headerdatainfo.push_back("RVmax_eig_zx");
        headerdatainfo.push_back("RVmax_eig_zy");
        headerdatainfo.push_back("RVmax_eig_zz");

        if (opt.iextrahalooutput) {
            headerdatainfo.push_back("Lx_200mean");
            headerdatainfo.push_back("Ly_200mean");
            headerdatainfo.push_back("Lz_200mean");
            headerdatainfo.push_back("Lx_200crit");
            headerdatainfo.push_back("Ly_200crit");
            headerdatainfo.push_back("Lz_200crit");
            headerdatainfo.push_back("Lx_BN98");
            headerdatainfo.push_back("Ly_BN98");
            headerdatainfo.push_back("Lz_BN98");
            if (opt.iInclusiveHalo>0) {
                headerdatainfo.push_back("Mass_200mean_excl");
                headerdatainfo.push_back("Mass_200crit_excl");
                headerdatainfo.push_back("Mass_BN98_excl");
                headerdatainfo.push_back("R_200mean_excl");
                headerdatainfo.push_back("R_200crit_excl");
                headerdatainfo.push_back("R_BN98_excl");
                headerdatainfo.push_back("Lx_200mean_excl");
                headerdatainfo.push_back("Ly_200mean_excl");
                headerdatainfo.push_back("Lz_200mean_excl");
                headerdatainfo.push_back("Lx_200crit_excl");
                headerdatainfo.push_back("Ly_200crit_excl");
                headerdatainfo.push_back("Lz_200crit_excl");
                headerdatainfo.push_back("Lx_BN98_excl");
                headerdatainfo.push_back("Ly_BN98_excl");
                headerdatainfo.push_back("Lz_BN98_excl");
            }
        }

#ifdef USEHDF
        sizeval=predtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif

#ifdef GASON
        headerdatainfo.push_back("n_gas");
#ifdef USEHDF
        predtypeinfo.push_back(PredType::STD_U64LE);
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
        headerdatainfo.push_back("M_gas");
        headerdatainfo.push_back("M_gas_Rvmax");
        headerdatainfo.push_back("M_gas_30kpc");
        //headerdatainfo.push_back("M_gas_50kpc");
        headerdatainfo.push_back("M_gas_500c");
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
    if (opt.iextragasoutput) {
        headerdatainfo.push_back("Mass_200mean_gas");
        headerdatainfo.push_back("Mass_200crit_gas");
        headerdatainfo.push_back("Mass_BN98_gas");
        headerdatainfo.push_back("Lx_200c_gas");
        headerdatainfo.push_back("Ly_200c_gas");
        headerdatainfo.push_back("Lz_200c_gas");
        headerdatainfo.push_back("Lx_200m_gas");
        headerdatainfo.push_back("Ly_200m_gas");
        headerdatainfo.push_back("Lz_200m_gas");
        headerdatainfo.push_back("Lx_BN98_gas");
        headerdatainfo.push_back("Ly_BN98_gas");
        headerdatainfo.push_back("Lz_BN98_gas");
        if (opt.iInclusiveHalo>0) {
            headerdatainfo.push_back("Mass_200mean_excl_gas");
            headerdatainfo.push_back("Mass_200crit_excl_gas");
            headerdatainfo.push_back("Mass_BN98_excl_gas");
            headerdatainfo.push_back("Lx_200c_excl_gas");
            headerdatainfo.push_back("Ly_200c_excl_gas");
            headerdatainfo.push_back("Lz_200c_excl_gas");
            headerdatainfo.push_back("Lx_200m_excl_gas");
            headerdatainfo.push_back("Ly_200m_excl_gas");
            headerdatainfo.push_back("Lz_200m_excl_gas");
            headerdatainfo.push_back("Lx_BN98_excl_gas");
            headerdatainfo.push_back("Ly_BN98_excl_gas");
            headerdatainfo.push_back("Lz_BN98_excl_gas");
        }
    }
#ifdef USEHDF
        sizeval=predtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif

#ifdef STARON
        headerdatainfo.push_back("n_star");
#ifdef USEHDF
        predtypeinfo.push_back(PredType::STD_U64LE);
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
        headerdatainfo.push_back("M_star");
        headerdatainfo.push_back("M_star_Rvmax");
        headerdatainfo.push_back("M_star_30kpc");
        //headerdatainfo.push_back("M_star_50kpc");
        headerdatainfo.push_back("M_star_500c");
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
        if (opt.iextrastaroutput) {
            headerdatainfo.push_back("Mass_200mean_star");
            headerdatainfo.push_back("Mass_200crit_star");
            headerdatainfo.push_back("Mass_BN98_star");
            headerdatainfo.push_back("Lx_200c_star");
            headerdatainfo.push_back("Ly_200c_star");
            headerdatainfo.push_back("Lz_200c_star");
            headerdatainfo.push_back("Lx_200m_star");
            headerdatainfo.push_back("Ly_200m_star");
            headerdatainfo.push_back("Lz_200m_star");
            headerdatainfo.push_back("Lx_BN98_star");
            headerdatainfo.push_back("Ly_BN98_star");
            headerdatainfo.push_back("Lz_BN98_star");
            if (opt.iInclusiveHalo>0) {
                headerdatainfo.push_back("Mass_200mean_excl_star");
                headerdatainfo.push_back("Mass_200crit_excl_star");
                headerdatainfo.push_back("Mass_BN98_excl_star");
                headerdatainfo.push_back("Lx_200c_excl_star");
                headerdatainfo.push_back("Ly_200c_excl_star");
                headerdatainfo.push_back("Lz_200c_excl_star");
                headerdatainfo.push_back("Lx_200m_excl_star");
                headerdatainfo.push_back("Ly_200m_excl_star");
                headerdatainfo.push_back("Lz_200m_excl_star");
                headerdatainfo.push_back("Lx_BN98_excl_star");
                headerdatainfo.push_back("Ly_BN98_excl_star");
                headerdatainfo.push_back("Lz_BN98_excl_star");
            }
        }
#ifdef USEHDF
        sizeval=predtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif

#ifdef BHON
        headerdatainfo.push_back("n_bh");
#ifdef USEHDF
        predtypeinfo.push_back(PredType::STD_U64LE);
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
        headerdatainfo.push_back("M_bh");
#ifdef USEHDF
        sizeval=predtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif


#ifdef HIGHRES
        headerdatainfo.push_back("n_interloper");
#ifdef USEHDF
        predtypeinfo.push_back(PredType::STD_U64LE);
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
        headerdatainfo.push_back("M_interloper");
#ifdef USEHDF
        sizeval=predtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif
        //additional information about a halo
        if (opt.iextrahalooutput) {

        }

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
    ///points to the the head pfof address of the group and parent
    Particle **Phead;
    Int_t **gidhead;
    ///parent pointers point to the address of the parents gidhead and Phead
    Particle **Pparenthead;
    Int_t **gidparenthead;
    ///add uber parent pointer (that is pointer to field halo)
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
    ///just allocate memory
    void Allocate(Int_t numgroups){
        nsinlevel=numgroups;
        Phead=new Particle*[numgroups+1];
        gidhead=new Int_t*[numgroups+1];
        stypeinlevel=new Int_t[numgroups+1];
        gidparenthead=new Int_t*[numgroups+1];
        giduberparenthead=new Int_t*[numgroups+1];
        nextlevel=NULL;
    }
    ///initialize
    void Initialize(){
        for (Int_t i=1;i<=nsinlevel;i++) {gidhead[i]=NULL;gidparenthead[i]=NULL;giduberparenthead[i]=NULL;}
    }
    ~StrucLevelData(){
        if (nextlevel!=NULL) delete nextlevel;
        nextlevel=NULL;
        if (nsinlevel>0) {
            delete[] Phead;
            delete[] gidhead;
            delete[] stypeinlevel;
            delete[] gidparenthead;
            delete[] giduberparenthead;
        }
    }
};

#if defined(USEHDF)||defined(USEADIOS)
///store the names of datasets in catalog output
struct DataGroupNames {
    ///store names of catalog group files
    vector<string> prop;
#ifdef USEHDF
    //store the data type
    vector<PredType> propdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospropdatatype;
#endif

    ///store names of catalog group files
    vector<string> group;
#ifdef USEHDF
    vector<PredType> groupdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiosgroupdatatype;
#endif

    ///store the names of catalog particle files
    vector<string> part;
#ifdef USEHDF
    vector<PredType> partdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospartdatatype;
#endif

    ///store the names of catalog particle files
    vector<string> types;
#ifdef USEHDF
    vector<PredType> typesdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiostypesdatatype;
#endif

    ///store the names of catalog particle files
    vector<string> hierarchy;
#ifdef USEHDF
    vector<PredType> hierarchydatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adioshierarchydatatype;
#endif

    ///store names of catalog group files
    vector<string> SO;
#ifdef USEHDF
    vector<PredType> SOdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> SOdatatype;
#endif

    DataGroupNames(){
        prop.push_back("File_id");
        prop.push_back("Num_of_files");
        prop.push_back("Num_of_groups");
        prop.push_back("Total_num_of_groups");
        prop.push_back("Cosmological_Sim");
        prop.push_back("Comoving_or_Physical");
        prop.push_back("Period");
        prop.push_back("Time");
        prop.push_back("Length_unit_to_kpc");
        prop.push_back("Velocity_to_kms");
        prop.push_back("Mass_unit_to_solarmass");
#ifdef USEHDF
        propdatatype.push_back(PredType::STD_I32LE);
        propdatatype.push_back(PredType::STD_I32LE);
        propdatatype.push_back(PredType::STD_U64LE);
        propdatatype.push_back(PredType::STD_U64LE);
        propdatatype.push_back(PredType::STD_U32LE);
        propdatatype.push_back(PredType::STD_U32LE);
        propdatatype.push_back(PredType::NATIVE_FLOAT);
        propdatatype.push_back(PredType::NATIVE_FLOAT);
        propdatatype.push_back(PredType::NATIVE_FLOAT);
        propdatatype.push_back(PredType::NATIVE_FLOAT);
        propdatatype.push_back(PredType::NATIVE_FLOAT);
#endif
#ifdef USEADIOS
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_real);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_real);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_real);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_real);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_real);
#endif

        group.push_back("File_id");
        group.push_back("Num_of_files");
        group.push_back("Num_of_groups");
        group.push_back("Total_num_of_groups");
        group.push_back("Group_Size");
        group.push_back("Offset");
        group.push_back("Offset_unbound");
#ifdef USEHDF
        groupdatatype.push_back(PredType::STD_I32LE);
        groupdatatype.push_back(PredType::STD_I32LE);
        groupdatatype.push_back(PredType::STD_U64LE);
        groupdatatype.push_back(PredType::STD_U64LE);
        groupdatatype.push_back(PredType::STD_U32LE);
        groupdatatype.push_back(PredType::STD_U64LE);
        groupdatatype.push_back(PredType::STD_U64LE);
#endif
#ifdef USEADIOS
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosgroupdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif

        part.push_back("File_id");
        part.push_back("Num_of_files");
        part.push_back("Num_of_particles_in_groups");
        part.push_back("Total_num_of_particles_in_all_groups");
        part.push_back("Particle_IDs");
#ifdef USEHDF
        partdatatype.push_back(PredType::STD_I32LE);
        partdatatype.push_back(PredType::STD_I32LE);
        partdatatype.push_back(PredType::STD_U64LE);
        partdatatype.push_back(PredType::STD_U64LE);
        partdatatype.push_back(PredType::STD_I64LE);
#endif
#ifdef USEADIOS
        adiospartdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiospartdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiospartdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospartdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospartdatatype.push_back(ADIOS_DATATYPES::adios_long);
#endif

        types.push_back("File_id");
        types.push_back("Num_of_files");
        types.push_back("Num_of_particles_in_groups");
        types.push_back("Total_num_of_particles_in_all_groups");
        types.push_back("Particle_types");
#ifdef USEHDF
        typesdatatype.push_back(PredType::STD_I32LE);
        typesdatatype.push_back(PredType::STD_I32LE);
        typesdatatype.push_back(PredType::STD_U64LE);
        typesdatatype.push_back(PredType::STD_U64LE);
        typesdatatype.push_back(PredType::STD_U16LE);
#endif
#ifdef USEADIOS
        adiostypesdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiostypesdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiostypesdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiostypesdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiostypesdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_short);
#endif

        hierarchy.push_back("File_id");
        hierarchy.push_back("Num_of_files");
        hierarchy.push_back("Num_of_groups");
        hierarchy.push_back("Total_num_of_groups");
        hierarchy.push_back("Number_of_substructures_in_halo");
        hierarchy.push_back("Parent_halo_ID");
#ifdef USEHDF
        hierarchydatatype.push_back(PredType::STD_I32LE);
        hierarchydatatype.push_back(PredType::STD_I32LE);
        hierarchydatatype.push_back(PredType::STD_U64LE);
        hierarchydatatype.push_back(PredType::STD_U64LE);
        hierarchydatatype.push_back(PredType::STD_U32LE);
        hierarchydatatype.push_back(PredType::STD_I64LE);
#endif
#ifdef USEADIOS
        adioshierarchydatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adioshierarchydatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adioshierarchydatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adioshierarchydatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adioshierarchydatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adioshierarchydatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
        SO.push_back("File_id");
        SO.push_back("Num_of_files");
        SO.push_back("Num_of_SO_regions");
        SO.push_back("Total_num_of_SO_regions");
        SO.push_back("Num_of_particles_in_SO_regions");
        SO.push_back("Total_num_of_particles_in_SO_regions");
        SO.push_back("SO_size");
        SO.push_back("Offset");
        SO.push_back("Particle_IDs");
#ifdef USEHDF
        SOdatatype.push_back(PredType::STD_I32LE);
        SOdatatype.push_back(PredType::STD_I32LE);
        SOdatatype.push_back(PredType::STD_U64LE);
        SOdatatype.push_back(PredType::STD_U64LE);
        SOdatatype.push_back(PredType::STD_U64LE);
        SOdatatype.push_back(PredType::STD_U64LE);
        SOdatatype.push_back(PredType::STD_U32LE);
        SOdatatype.push_back(PredType::STD_U64LE);
        SOdatatype.push_back(PredType::STD_I64LE);
#endif
#ifdef USEADIOS
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_long);
#endif
    }
};
#endif

///if using MPI API
#ifdef USEMPI
#include <mpi.h>
///Includes external global variables used for MPI version of code
#include "mpivar.h"
#endif

extern StrucLevelData *psldata;

#endif
