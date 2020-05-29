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
#include <set>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <unistd.h>

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
#include "ompvar.h"

///if using HDF API
#ifdef USEHDF
#include "hdf5.h"
#endif

///if using ADIOS API
#ifdef USEADIOS
#include "adios.h"
#endif

#include "git_revision.h"

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
#define SUBSTYPE 10
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
///baryon 6D FOF search
#define FOFBARYON6D 0
///baryon phase tensor search
#define FOFBARYONPHASETENSOR 1
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

///\defgroup CALCULATIONTYPES defining what is calculated
//@{
#define CALCAVERAGE 1
#define CALCTOTAL 2
#define CALCSTD 3
#define CALCMEDIAN 4
#define CALCMIN 5
#define CALCMAX 6
#define CALCLOGAVERAGE 7
#define CALCLOGSTD 8
#define CALCQUANTITYMASSWEIGHT 10
#define CALCAVERAGEMASSWEIGHT 11
#define CALCTOTALMASSWEIGHT 12
#define CALCSTDMASSWEIGHT 13
#define CALCMEDIANMASSWEIGHT 14
#define CALCMINMASSWEIGHT 15
#define CALCMAXMASSWEIGHT 16
#define CALCLOGAVERAGEMASSWEIGHT 17
#define CALCLOGSTDMASSWEIGHT 18
#define CALCQUANTITYAPERTURETOTAL -1
#define CALCQUANTITYAPERTUREAVERAGE -2
typedef double (*ExtraPropFunc)(double, double, double&);

//@}

/// \name For Unbinding
//@{

///number below which just use PP calculation for potential, which occurs roughly at when n~2*log(n) (from scaling of n^2 vs n ln(n) for PP vs tree and factor of 2 is
///for extra overhead in producing tree. For reasonable values of n (>100) this occurs at ~100. Here to account for extra memory need for tree, we use n=3*log(n) or 150
#define UNBINDNUM 150
#define POTPPCALCNUM 150
#define POTOMPCALCNUM 1000
///when unbinding check to see if system is bound and least bound particle is also bound
#define USYSANDPART 0
///when unbinding check to see if least bound particle is also bound
#define UPART 1
///use the bulk centre of mass velocity to define velocity reference frame when determining if particle bound
#define CMVELREF 0
///use the particle at potential minimum. Issues if too few particles used as particles will move in and out of deepest point of the potential well
#define POTREF 1
///use Centre-of-mass to caculate Properties
#define PROPREFCM 0
///use most bound particle to calculate properties
#define PROPREFMBP 1
///use minimum potential particle to calculat properties
#define PROPREFMINPOT 2

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

/// \defgroup PROPLIMS Particle limits for calculating properties
//@{
#define PROPNFWMINNUM 100
#define PROPCMMINNUM 10
#define PROPROTMINNUM 10
#define PROPMORPHMINNUM 10
//@}

///\name halo id modifers used with current snapshot value to make temporally unique halo identifiers
#ifdef LONGINT
#define HALOIDSNVAL 1000000000000L
#else
#define HALOIDSNVAL 1000000
#endif

///\defgroup radial profile parameters
//@{
#define PROFILERNORMPHYS 0
#define PROFILERNORMR200CRIT 1
#define PROFILERBINTYPELOG 0
//@}

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
    ///Maximum fraction of particles that can be considered unbound before group removed entirely
    Double_t maxunboundfracforiterativeunbind;
    ///Max allowed unbound fraction to speed up unbinding
    Double_t maxallowedunboundfrac;

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
        TreeThetaOpen=0.5;
        eps=0.0;
        Npotref=20;
        fracpotref=1.0;
        maxunbindfrac=0.5;
        maxunboundfracforiterativeunbind=0.95;
        maxallowedunboundfrac=0.025;
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

    ///\name git related info
    //@{
    string git_sha1;
    //@}
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
    /// calculate and output extra bh fields
    int iextrabhoutput;
    /// calculate and output extra interloper fields
    int iextrainterloperoutput;
    /// calculate subind like properties
    int isubfindproperties;
    ///for output, produce subfind like format
    int isubfindoutput;
    ///flag indicating that VR is running on the fly
    bool iontheflyfinding;

    ///disable particle id related output like fof.grp or catalog_group data. Useful if just want halo properties
    ///and not interested in tracking. Code writes halo properties catalog and exits.
    int inoidoutput;
    ///return propery data in in comoving little h units instead of standard physical units
    int icomoveunit;
    /// input is a cosmological simulation so can use box sizes, cosmological parameters, etc to set scales
    int icosmologicalin;
    /// input buffer size when reading data
    long long inputbufsize;
    /// mpi paritcle buffer size when sending input particle information
    long long mpiparticletotbufsize,mpiparticlebufsize;
    /// mpi factor by which to multiple the memory allocated, ie: buffer region
    /// to reduce likelihood of having to expand/allocate new memory
    Double_t mpipartfac;
    /// if using parallel output, number of mpi threads to group together
    int mpinprocswritesize;

    /// run FOF using OpenMP
    int iopenmpfof;
    /// size of openmp FOF region
    int openmpfofsize;

    ///\name length,m,v,grav conversion units
    //@{
    Double_t lengthinputconversion, massinputconversion, energyinputconversion, internalenergyinputconversion, velocityinputconversion;
    Double_t SFRinputconversion, metallicityinputconversion, stellarageinputconversion;
    int istellaragescalefactor, isfrisssfr;
    Double_t G;
    Double_t lengthtokpc, velocitytokms, masstosolarmass, energyperunitmass, timetoseconds;
    Double_t SFRtosolarmassperyear, stellaragetoyrs, metallicitytosolar;
    //@}
    ///period (comove)
    Double_t p;
    ///\name scale factor, Hubunit, h, cosmology, virial density. These are used if linking lengths are scaled or trying to define virlevel using the cosmology
    //@{
    Double_t a,H,h;
    Double_t Omega_m, Omega_b, Omega_cdm, Omega_Lambda, Omega_k, Omega_r, Omega_nu, Omega_de, w_de;
    Double_t rhocrit, rhobg, virlevel, virBN98;
    int comove;
    /// to store the internal code unit to kpc and the distance^2 of 30 kpc, and 50 kpc
    Double_t lengthtokpc30pow2, lengthtokpc50pow2;
    //@}

    ///to store number of each particle types so that if only searching one particle type, assumes order of gas, dark (halo, disk,bulge), star, special or sink for pfof tipsy style output
    Int_t numpart[NPARTTYPES];

    ///\name parameters that control the local and average volumes used to calculate the local velocity density and the mean field, also the size of the leafnode in the kd-tree used when searching the tree for fof neighbours
    //@{
    int iLocalVelDenApproxCalcFlag;
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
    /// FOF search for baryons
    int ifofbaryonsearch;
    ///flag indicating if move to CM frame for substructure search
    int icmrefadjust;
    /// flag indicating if CM is interated shrinking spheres
    int iIterateCM;
    /// flag to sort output particle lists by binding energy (or potential if not on)
    int iSortByBindingEnergy;
    /// what reference position to use when calculating Properties
    int iPropertyReferencePosition;
    /// what particle type is used to define reference position
    int ParticleTypeForRefenceFrame;


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
    Double_t ellhalo3dxfac;
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
    ///factor by which a substructure s must be closer than in phase-space to merger with another substructure in sigma units
    Double_t coresubmergemindist;
    ///whether substructure phase-space distance merge check is applied to background host halo as well.
    int icoresubmergewithbg;
    ///fraction of size a substructure must be of host to be considered a spurious dynamical substructure
    Double_t minfracsubsizeforremoval;
    ///Maximum allowed mean local velocity density ratio above which structure is considered highly unrelaxed.
    Double_t maxmeanlocalvelratio;
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
    /// input contains dm particles
    int iusedmparticles;
    /// input contains hydro/gas particles
    int iusegasparticles;
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
    ///minimum enclosed mass on which to base SO calculations, <1
    Double_t SphericalOverdensityMinHaloFac;
    ///if want to the particle IDs that are within the SO overdensity of a halo
    int iSphericalOverdensityPartList;
    /// if want to include more than just field objects (halos) in full SO calculations
    int SphericalOverdensitySeachMaxStructLevel;
    /// flag to store whether SO calculations need extra properties
    bool iSphericalOverdensityExtraFieldCalculations;
    /// \name Extra variables to store information useful in zoom simluations
    //@{
    /// store the lowest dark matter particle mass
    Double_t zoomlowmassdm;
    //@}

    ///\name extra runtime flags
    //@{
    ///scale lengths. Useful if searching single halo system and which to automatically scale linking lengths
    int iScaleLengths;

    /// \name Swift/Metis related quantitites
    //@{
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

    /// \name options related to calculation of aperture/profile
    //@{
    int iaperturecalc;
    int aperturenum,apertureprojnum;
    vector<Double_t> aperture_values_kpc;
    vector<string> aperture_names_kpc;
    vector<Double_t> aperture_proj_values_kpc;
    vector<string> aperture_proj_names_kpc;
    int iprofilecalc, iprofilenorm, iprofilebintype;
    int profilenbins;
    int iprofilecumulative;
    string profileradnormstring;
    vector<Double_t> profile_bin_edges;
    Int_t profileminsize, profileminFOFsize;
    //@}

    /// \name options related to calculation of arbitrary overdensities masses, radii, angular momentum
    //@{
    int SOnum;
    vector<Double_t> SOthresholds_values_crit;
    vector<string> SOthresholds_names_crit;
    //@}

    /// \name options related to calculating star forming gas quantities
    //@{
    Double_t gas_sfr_threshold;
    //@}

    /// \name options related to calculating detailed hydro/star/bh properties related to chemistry/feedbac, etc
    //@{
    ///stores the name of the field
    vector<string> gas_internalprop_names;
    vector<string> star_internalprop_names;
    vector<string> bh_internalprop_names;

    ///can also store the dimensional index of the field, useful when single data
    ///set contains many related but different properties
    vector<unsigned int> gas_internalprop_index;
    vector<unsigned int> star_internalprop_index;
    vector<unsigned int> bh_internalprop_index;
    ///stores what is calculated
    ///(1 is mass weighted average, 2 mass weighted total, etc
    vector<int> gas_internalprop_function;
    vector<int> star_internalprop_function;
    vector<int> bh_internalprop_function;

    vector<string> gas_chem_names;
    vector<string> star_chem_names;
    vector<string> bh_chem_names;
    vector<unsigned int> gas_chem_index;
    vector<unsigned int> star_chem_index;
    vector<unsigned int> bh_chem_index;
    vector<int> gas_chem_function;
    vector<int> star_chem_function;
    vector<int> bh_chem_function;

    vector<string> gas_chemproduction_names;
    vector<string> star_chemproduction_names;
    vector<string> bh_chemproduction_names;
    vector<unsigned int> gas_chemproduction_index;
    vector<unsigned int> star_chemproduction_index;
    vector<unsigned int> bh_chemproduction_index;
    vector<int> gas_chemproduction_function;
    vector<int> star_chemproduction_function;
    vector<int> bh_chemproduction_function;

    vector<string> extra_dm_internalprop_names;
    vector<unsigned int> extra_dm_internalprop_index;
    vector<int> extra_dm_internalprop_function;

    ///store the output field name
    vector<string> gas_internalprop_output_names;
    vector<string> star_internalprop_output_names;
    vector<string> bh_internalprop_output_names;
    vector<string> gas_chem_output_names;
    vector<string> star_chem_output_names;
    vector<string> bh_chem_output_names;
    vector<string> gas_chemproduction_output_names;
    vector<string> star_chemproduction_output_names;
    vector<string> bh_chemproduction_output_names;
    vector<string> extra_dm_internalprop_output_names;

    ///store conversion factor from input unit to output unit
    vector<float> gas_internalprop_input_output_unit_conversion_factors;
    vector<float> star_internalprop_input_output_unit_conversion_factors;
    vector<float> bh_internalprop_input_output_unit_conversion_factors;
    vector<float> gas_chem_input_output_unit_conversion_factors;
    vector<float> star_chem_input_output_unit_conversion_factors;
    vector<float> bh_chem_input_output_unit_conversion_factors;
    vector<float> gas_chemproduction_input_output_unit_conversion_factors;
    vector<float> star_chemproduction_input_output_unit_conversion_factors;
    vector<float> bh_chemproduction_input_output_unit_conversion_factors;
    vector<float> extra_dm_internalprop_input_output_unit_conversion_factors;

    ///store output units
    vector<string> gas_internalprop_output_units;
    vector<string> star_internalprop_output_units;
    vector<string> bh_internalprop_output_units;
    vector<string> gas_chem_output_units;
    vector<string> star_chem_output_units;
    vector<string> bh_chem_output_units;
    vector<string> gas_chemproduction_output_units;
    vector<string> star_chemproduction_output_units;
    vector<string> bh_chemproduction_output_units;
    vector<string> extra_dm_internalprop_output_units;

    ///some calculations are multistage and must be paired with
    ///another calculation in the list. This is true of standard deviations
    vector<int> gas_internalprop_index_paired_calc;
    vector<int> star_internalprop_index_paired_calc;
    vector<int> bh_internalprop_index_paired_calc;
    vector<int> gas_chem_index_paired_calc;
    vector<int> star_chem_index_paired_calc;
    vector<int> bh_chem_index_paired_calc;
    vector<int> gas_chemproduction_index_paired_calc;
    vector<int> star_chemproduction_index_paired_calc;
    vector<int> bh_chemproduction_index_paired_calc;
    vector<int> extra_dm_internalprop_index_paired_calc;

    ///whether some calculations are for extra properties are aperture calculations
    bool gas_extraprop_aperture_calc;
    bool star_extraprop_aperture_calc;
    bool bh_extraprop_aperture_calc;
    bool extra_dm_extraprop_aperture_calc;

    ///easier to store information separately internally
    vector<string> gas_internalprop_names_aperture;
    vector<string> gas_chem_names_aperture;
    vector<string> gas_chemproduction_names_aperture;
    vector<string> star_internalprop_names_aperture;
    vector<string> star_chem_names_aperture;
    vector<string> star_chemproduction_names_aperture;
    vector<string> bh_internalprop_names_aperture;
    vector<string> bh_chem_names_aperture;
    vector<string> bh_chemproduction_names_aperture;
    vector<string> extra_dm_internalprop_names_aperture;

    vector<unsigned int> gas_internalprop_index_aperture;
    vector<unsigned int> gas_chem_index_aperture;
    vector<unsigned int> gas_chemproduction_index_aperture;
    vector<unsigned int> star_internalprop_index_aperture;
    vector<unsigned int> star_chem_index_aperture;
    vector<unsigned int> star_chemproduction_index_aperture;
    vector<unsigned int> bh_internalprop_index_aperture;
    vector<unsigned int> bh_chem_index_aperture;
    vector<unsigned int> bh_chemproduction_index_aperture;
    vector<unsigned int> extra_dm_internalprop_index_aperture;

    vector<int> gas_internalprop_function_aperture;
    vector<int> gas_chem_function_aperture;
    vector<int> gas_chemproduction_function_aperture;
    vector<int> star_internalprop_function_aperture;
    vector<int> star_chem_function_aperture;
    vector<int> star_chemproduction_function_aperture;
    vector<int> bh_internalprop_function_aperture;
    vector<int> bh_chem_function_aperture;
    vector<int> bh_chemproduction_function_aperture;
    vector<int> extra_dm_internalprop_function_aperture;

    vector<string> gas_internalprop_output_units_aperture;
    vector<string> star_internalprop_output_units_aperture;
    vector<string> bh_internalprop_output_units_aperture;
    vector<string> gas_chem_output_units_aperture;
    vector<string> star_chem_output_units_aperture;
    vector<string> bh_chem_output_units_aperture;
    vector<string> gas_chemproduction_output_units_aperture;
    vector<string> star_chemproduction_output_units_aperture;
    vector<string> bh_chemproduction_output_units_aperture;
    vector<string> extra_dm_internalprop_output_units_aperture;

    vector<float> gas_internalprop_input_output_unit_conversion_factors_aperture;
    vector<float> star_internalprop_input_output_unit_conversion_factors_aperture;
    vector<float> bh_internalprop_input_output_unit_conversion_factors_aperture;
    vector<float> gas_chem_input_output_unit_conversion_factors_aperture;
    vector<float> star_chem_input_output_unit_conversion_factors_aperture;
    vector<float> bh_chem_input_output_unit_conversion_factors_aperture;
    vector<float> gas_chemproduction_input_output_unit_conversion_factors_aperture;
    vector<float> star_chemproduction_input_output_unit_conversion_factors_aperture;
    vector<float> bh_chemproduction_input_output_unit_conversion_factors_aperture;
    vector<float> extra_dm_internalprop_input_output_unit_conversion_factors_aperture;

    vector<string> gas_internalprop_output_names_aperture;
    vector<string> gas_chem_output_names_aperture;
    vector<string> gas_chemproduction_output_names_aperture;
    vector<string> star_internalprop_output_names_aperture;
    vector<string> star_chem_output_names_aperture;
    vector<string> star_chemproduction_output_names_aperture;
    vector<string> bh_internalprop_output_names_aperture;
    vector<string> bh_chem_output_names_aperture;
    vector<string> bh_chemproduction_output_names_aperture;
    vector<string> extra_dm_internalprop_output_names_aperture;

    //to store the unique names that are going to be loaded from the input
    vector<string> gas_internalprop_unique_input_names;
    vector<string> gas_chem_unique_input_names;
    vector<string> gas_chemproduction_unique_input_names;
    vector<string> star_internalprop_unique_input_names;
    vector<string> star_chem_unique_input_names;
    vector<string> star_chemproduction_unique_input_names;
    vector<string> bh_internalprop_unique_input_names;
    vector<string> bh_chem_unique_input_names;
    vector<string> bh_chemproduction_unique_input_names;
    vector<string> extra_dm_internalprop_unique_input_names;

    vector<unsigned short> gas_internalprop_unique_input_indexlist;
    vector<unsigned short> gas_chem_unique_input_indexlist;
    vector<unsigned short> gas_chemproduction_unique_input_indexlist;
    vector<unsigned short> star_internalprop_unique_input_indexlist;
    vector<unsigned short> star_chem_unique_input_indexlist;
    vector<unsigned short> star_chemproduction_unique_input_indexlist;
    vector<unsigned short> bh_internalprop_unique_input_indexlist;
    vector<unsigned short> bh_chem_unique_input_indexlist;
    vector<unsigned short> bh_chemproduction_unique_input_indexlist;
    vector<unsigned short> extra_dm_internalprop_unique_input_indexlist;

    //@}

    /// \name memory related info
    //@{
    unsigned long long memuse_peak;
    unsigned long long memuse_ave;
    int memuse_nsamples;
    bool memuse_log;
    //@}

    //silly flag to store whether input has little h's in it.
    bool inputcontainslittleh;

    Options()
    {
        lengthinputconversion = 1.0;
        massinputconversion = 1.0;
        velocityinputconversion = 1.0;
        SFRinputconversion = 1.0;
        metallicityinputconversion = 1.0;
        energyinputconversion = 1.0;
        stellarageinputconversion =1.0;
        istellaragescalefactor = 1;
        isfrisssfr = 0;

        G = 1.0;
        p = 0.0;

        a = 1.0;
        H = 0.;
        h = 1.0;
        Omega_m = 1.0;
        Omega_Lambda = 0.0;
        Omega_b = 0.0;
        Omega_cdm = Omega_m;
        Omega_k = 0;
        Omega_r = 0.0;
        Omega_nu = 0.0;
        Omega_de = 0.0;
        w_de = -1.0;
        rhobg = 1.0;
        virlevel = -1;
        comove=0;
        H=100.0;//value of Hubble flow in h 1 km/s/Mpc
        MassValue=1.0;

        inputtype=IOGADGET;

        num_files=1;
        nsnapread=1;

        fname=outname=smname=pname=gname=outname=NULL;

        Bsize=32;
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
        iIterateCM = 1;
        iLocalVelDenApproxCalcFlag = 2 ;

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
        ellhalo3dxfac=-1.0;

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
        iSortByBindingEnergy=1;
        iPropertyReferencePosition=PROPREFCM;
        ParticleTypeForRefenceFrame=-1;

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
        coresubmergemindist=0.0;
        icoresubmergewithbg=0;
        minfracsubsizeforremoval=0.75;
        maxmeanlocalvelratio=0.5;

        iverbose=0;
        iwritefof=0;
        iseparatefiles=0;
        ibinaryout=0;
        iextendedoutput=0;
        isubfindoutput=0;
        inoidoutput=0;
        icomoveunit=0;
        icosmologicalin=1;

        iextrahalooutput=0;
        iextragasoutput=0;
        iextrastaroutput=0;
        iextrainterloperoutput=0;
        isubfindproperties=0;

        iusedmparticles=1;
        iusegasparticles=1;
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

        inputbufsize=1000000;

        mpiparticletotbufsize=-1;
        mpiparticlebufsize=-1;
        mpinprocswritesize=1;

        lengthtokpc=-1.0;
        velocitytokms=-1.0;
        masstosolarmass=-1.0;
#if defined(GASON) || defined(STARON) || defined(BHON)
        SFRtosolarmassperyear=-1.0;
        stellaragetoyrs=-1.0;
        metallicitytosolar=-1.0;
#endif


        lengthtokpc30pow2=30.0*30.0;
        lengthtokpc50pow2=50.0*50.0;

        SphericalOverdensitySeachFac=2.5;
        SphericalOverdensityMinHaloFac=0.05;
        iSphericalOverdensityPartList=0;
        SphericalOverdensitySeachMaxStructLevel = HALOSTYPE;
        iSphericalOverdensityExtraFieldCalculations = false;

        mpipartfac=0.1;
#if USEHDF
        ihdfnameconvention=-1;
#endif
        iaperturecalc=0;
        aperturenum=0;
        apertureprojnum=0;
        SOnum=0;

        iprofilecalc=0;
        iprofilenorm=PROFILERNORMR200CRIT;
        iprofilebintype=PROFILERBINTYPELOG;
        iprofilecumulative=0;
        profilenbins=0;
        profileminsize = profileminFOFsize = 0;
#ifdef USEOPENMP
        iopenmpfof = 1;
        openmpfofsize = ompfofsearchnum;
#endif

        iontheflyfinding = false;

        memuse_peak = 0;
        memuse_ave = 0;
        memuse_nsamples = 0;
        memuse_log = false;

        inputcontainslittleh = true;

    }
    Options(Options &opt) = default;
    Options& operator=(const Options&) = default;
    Options& operator=(Options&&) = default;
};

struct ConfigInfo{
    //list the name of the info
    vector<string> nameinfo;
    vector<string> datainfo;
    vector<string> datatype;

    string python_type_string(bool &x){return string("bool");}
    string python_type_string(int &x){return string("int32");}
    string python_type_string(unsigned int &x){return string("uint32");}
    string python_type_string(long long &x){return string("int64");}
    string python_type_string(unsigned long long &x){return string("uint64");}
    string python_type_string(float &x){return string("float32");}
    string python_type_string(double &x){return string("float64");}
    string python_type_string(string &x){return string("str");}

    void AddEntry(string entryname){
        nameinfo.push_back(entryname);
        datainfo.push_back("");
        datatype.push_back("");
    }
    template<typename T> void AddEntry(string entryname, T entry){
        nameinfo.push_back(entryname);
        datainfo.push_back(to_string(entry));
        datatype.push_back(python_type_string(entry));
    }
    template<typename T> void AddEntry(string entryname, vector<T> entries){
        if (entries.size() == 0) return;
        T val = entries[0];
        nameinfo.push_back(entryname);
        string datastring=string("");
        for (auto &x:entries) {datastring+=to_string(x);datastring+=string(",");}
        datainfo.push_back(datastring);
        datatype.push_back(python_type_string(val));
    }
    template<typename T> void AddEntry(string entryname, vector<T> entries1, vector<T> entries2){
        if (entries1.size() + entries2.size() == 0) return;
        T val = entries1[0];
        nameinfo.push_back(entryname);
        string datastring=string("");
        for (auto &x:entries1) {datastring+=to_string(x);datastring+=string(",");}
        for (auto &x:entries2) {datastring+=to_string(x);datastring+=string(",");}
        datainfo.push_back(datastring);
        datatype.push_back(python_type_string(val));
    }
    void AddEntry(string entryname, string entry){
        nameinfo.push_back(entryname);
        datainfo.push_back(entry);
        datatype.push_back(python_type_string(entry));
    }
    void AddEntry(string entryname, vector<string> entries){
        if (entries.size() == 0) return;
        string val = entries[0];
        nameinfo.push_back(entryname);
        string datastring=string("");
        for (auto &x:entries) {datastring+=x;datastring+=string(",");}
        datainfo.push_back(datastring);
        datatype.push_back(python_type_string(val));
    }
    void AddEntry(string entryname, vector<string> entries1, vector<string> entries2){
        if (entries1.size() + entries2.size() == 0) return;
        string val = entries1[0];
        nameinfo.push_back(entryname);
        string datastring=string("");
        for (auto &x:entries1) {datastring+=x;datastring+=string(",");}
        for (auto &x:entries2) {datastring+=x;datastring+=string(",");}
        datainfo.push_back(datastring);
        datatype.push_back(python_type_string(val));
    }

    ConfigInfo(Options &opt);
};

struct SimInfo{
    //list the name of the info
    vector<string> nameinfo;
    vector<string> datainfo;
    vector<string> datatype;

    string python_type_string(int &x){return string("int32");}
    string python_type_string(unsigned int &x){return string("uint32");}
    string python_type_string(long &x){return string("int64");}
    string python_type_string(unsigned long &x){return string("uint64");}
    string python_type_string(float &x){return string("float32");}
    string python_type_string(double &x){return string("float64");}

    SimInfo(Options &opt){
        //if compiler is super old and does not have at least std 11 implementation to_string does not exist
#ifndef OLDCCOMPILER
        nameinfo.push_back("Cosmological_Sim");
        datainfo.push_back(to_string(opt.icosmologicalin));
        datatype.push_back(python_type_string(opt.icosmologicalin));
        if (opt.icosmologicalin) {
            nameinfo.push_back("ScaleFactor");
            datainfo.push_back(to_string(opt.a));
            datatype.push_back(python_type_string(opt.a));
            nameinfo.push_back("h_val");
            datainfo.push_back(to_string(opt.h));
            datatype.push_back(python_type_string(opt.h));
            nameinfo.push_back("Omega_m");
            datainfo.push_back(to_string(opt.Omega_m));
            datatype.push_back(python_type_string(opt.Omega_m));
            nameinfo.push_back("Omega_Lambda");
            datainfo.push_back(to_string(opt.Omega_Lambda));
            datatype.push_back(python_type_string(opt.Omega_Lambda));
            nameinfo.push_back("Omega_cdm");
            datainfo.push_back(to_string(opt.Omega_cdm));
            datatype.push_back(python_type_string(opt.Omega_cdm));
            nameinfo.push_back("Omega_b");
            datainfo.push_back(to_string(opt.Omega_b));
            datatype.push_back(python_type_string(opt.Omega_b));
            nameinfo.push_back("w_of_DE");
            datainfo.push_back(to_string(opt.w_de));
            datatype.push_back(python_type_string(opt.w_de));
            nameinfo.push_back("Period");
            datainfo.push_back(to_string(opt.p));
            datatype.push_back(python_type_string(opt.p));
            nameinfo.push_back("Hubble_unit");
            datainfo.push_back(to_string(opt.H));
            datatype.push_back(python_type_string(opt.H));
        }
        else{
            nameinfo.push_back("Time");
            datainfo.push_back(to_string(opt.a));
            datatype.push_back(python_type_string(opt.a));
            nameinfo.push_back("Period");
            datainfo.push_back(to_string(opt.p));
            datatype.push_back(python_type_string(opt.p));
        }

        //units
        nameinfo.push_back("Length_unit");
        datainfo.push_back(to_string(opt.lengthinputconversion));
        datatype.push_back(python_type_string(opt.lengthinputconversion));
        nameinfo.push_back("Velocity_unit");
        datainfo.push_back(to_string(opt.velocityinputconversion));
        datatype.push_back(python_type_string(opt.velocityinputconversion));
        nameinfo.push_back("Mass_unit");
        datainfo.push_back(to_string(opt.massinputconversion));
        datatype.push_back(python_type_string(opt.massinputconversion));
        nameinfo.push_back("Gravity");
        datainfo.push_back(to_string(opt.G));
        datatype.push_back(python_type_string(opt.G));
#ifdef NOMASS
        nameinfo.push_back("Mass_value");
        datainfo.push_back(to_string(opt.MassValue));
        datatype.push_back(python_type_string(opt.MassValue));
#endif

#endif
    }
};


struct UnitInfo{
    //list the name of the info
    vector<string> nameinfo;
    vector<string> datainfo;
    vector<string> datatype;

    string python_type_string(int &x){return string("int32");}
    string python_type_string(unsigned int &x){return string("uint32");}
    string python_type_string(long &x){return string("int64");}
    string python_type_string(unsigned long &x){return string("uint64");}
    string python_type_string(float &x){return string("float32");}
    string python_type_string(double &x){return string("float64");}

    UnitInfo(Options &opt){
        //if compiler is super old and does not have at least std 11 implementation to_string does not exist
#ifndef OLDCCOMPILER
        nameinfo.push_back("Cosmological_Sim");
        datainfo.push_back(to_string(opt.icosmologicalin));
        datatype.push_back(python_type_string(opt.icosmologicalin));
        nameinfo.push_back("Comoving_or_Physical");
        datainfo.push_back(to_string(opt.icomoveunit));
        datatype.push_back(python_type_string(opt.icomoveunit));
        //units
        nameinfo.push_back("Length_unit_to_kpc");
        datainfo.push_back(to_string(opt.lengthtokpc));
        datatype.push_back(python_type_string(opt.lengthtokpc));
        nameinfo.push_back("Velocity_unit_to_kms");
        datainfo.push_back(to_string(opt.velocitytokms));
        datatype.push_back(python_type_string(opt.velocitytokms));
        nameinfo.push_back("Mass_unit_to_solarmass");
        datainfo.push_back(to_string(opt.masstosolarmass));
        datatype.push_back(python_type_string(opt.masstosolarmass));
#if defined(GASON) || defined(STARON) || defined(BHON)
        nameinfo.push_back("Metallicity_unit_to_solar");
        datainfo.push_back(to_string(opt.metallicitytosolar));
        datatype.push_back(python_type_string(opt.metallicitytosolar));
        nameinfo.push_back("SFR_unit_to_solarmassperyear");
        datainfo.push_back(to_string(opt.SFRtosolarmassperyear));
        datatype.push_back(python_type_string(opt.SFRtosolarmassperyear));
        nameinfo.push_back("Stellar_age_unit_to_yr");
        datainfo.push_back(to_string(opt.stellaragetoyrs));
        datatype.push_back(python_type_string(opt.stellaragetoyrs));
#endif
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
        if (nparts>0)delete[] nindex;
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
    Int_t gNFOF,gN6DFOF;
    ///centre of mass
    Coordinate gcm, gcmvel;
    ///Position of most bound particle, and also of particle with min potential
    Coordinate gposmbp, gvelmbp, gposminpot, gvelminpot;
    ///\name physical properties regarding mass, size
    //@{
    Double_t gmass,gsize,gMvir,gRvir,gRcm,gRmbp,gRminpot,gmaxvel,gRmaxvel,gMmaxvel,gRhalfmass,gMassTwiceRhalfmass;
    Double_t gM200c,gR200c,gM200m,gR200m,gMFOF,gM6DFOF,gM500c,gR500c,gMBN98,gRBN98;
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
    ///Keep track of position of least unbound particle and most bound particle pid and minimum potential
    Int_t iunbound,ibound, iminpot;
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

    ///\name radial profiles
    //@{
    vector<unsigned int> aperture_npart;
    vector<float> aperture_mass;
    vector<float> aperture_veldisp;
    vector<float> aperture_vrdisp;
    vector<float> aperture_rhalfmass;
    vector<Coordinate> aperture_mass_proj;
    vector<Coordinate> aperture_rhalfmass_proj;
    vector<Coordinate> aperture_L;
    vector<unsigned int> profile_npart;
    vector<unsigned int> profile_npart_inclusive;
    vector<float> profile_mass;
    vector<float> profile_mass_inclusive;
    vector<Coordinate> profile_L;
    #if defined(GASON) || defined(STARON) || defined(BHON)
    vector<unsigned int> aperture_npart_dm;
    vector<float> aperture_mass_dm;
    vector<float> aperture_veldisp_dm;
    vector<float> aperture_vrdisp_dm;
    vector<float> aperture_rhalfmass_dm;
    #endif
    //@}

    vector<Double_t> SO_mass, SO_radius;
    vector<Coordinate> SO_angularmomentum;

#ifdef GASON
    ///\name gas specific quantities
    //@{
    ///number of particles
    int n_gas;
    ///mass
    Double_t M_gas, M_gas_rvmax, M_gas_30kpc, M_gas_50kpc, M_gas_500c;
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
    Double_t MassTwiceRhalfmass_gas, Rhalfmass_gas, q_gas, s_gas;
    Matrix eigvec_gas;
    ///mass weighted sum of temperature, metallicty, star formation rate
    Double_t Temp_gas, Z_gas, SFR_gas;
    ///mean temperature,metallicty,star formation rate
    Double_t Temp_mean_gas, Z_mean_gas, SFR_mean_gas;
    ///physical properties for dynamical state
    Double_t Efrac_gas, Pot_gas, T_gas;
    //@}

    ///\name gas radial profiles
    //@{
    vector<unsigned int> aperture_npart_gas;
    vector<float> aperture_mass_gas;
    vector<float> aperture_veldisp_gas;
    vector<float> aperture_vrdisp_gas;
    vector<float> aperture_SFR_gas;
    vector<float> aperture_Z_gas;
    vector<float> aperture_rhalfmass_gas;
    vector<Coordinate> aperture_L_gas;
    vector<Coordinate> aperture_mass_proj_gas;
    vector<Coordinate> aperture_rhalfmass_proj_gas;
    vector<Coordinate> aperture_SFR_proj_gas;
    vector<Coordinate> aperture_Z_proj_gas;
    vector<unsigned int> profile_npart_gas;
    vector<unsigned int> profile_npart_inclusive_gas;
    vector<float> profile_mass_gas;
    vector<float> profile_mass_inclusive_gas;
    vector<Coordinate> profile_L_gas;
    //@}

    vector<Double_t> SO_mass_gas;
    vector<Coordinate> SO_angularmomentum_gas;
#ifdef STARON
    ///\name star forming gas specific quantities
    //@{
    ///number of particles
    int n_gas_sf;
    ///mass
    Double_t M_gas_sf, M_gas_sf_rvmax,M_gas_sf_30kpc,M_gas_sf_50kpc, M_gas_sf_500c;
    ///mass in spherical overdensities
    Double_t M_200crit_gas_sf, M_200mean_gas_sf, M_BN98_gas_sf;
    ///mass in spherical overdensities inclusive of all masses
    Double_t M_200crit_excl_gas_sf, M_200mean_excl_gas_sf, M_BN98_excl_gas_sf;
    ///velocity/angular momentum info
    Double_t Krot_gas_sf;
    Coordinate L_gas_sf;
    ///physical properties for angular momentum (can be inclusive or exclusive )
    Coordinate L_200crit_gas_sf, L_200mean_gas_sf, L_BN98_gas_sf;
    ///physical properties for angular momentum exclusiveto object
    Coordinate L_200crit_excl_gas_sf, L_200mean_excl_gas_sf, L_BN98_excl_gas_sf;
    //dispersion
    Double_t sigV_gas_sf;
    ///morphology
    Double_t MassTwiceRhalfmass_gas_sf, Rhalfmass_gas_sf, q_gas_sf, s_gas_sf;
    ///mass weighted sum of temperature, metallicty, star formation rate
    Double_t Temp_gas_sf, Z_gas_sf, SFR_gas_sf;
    ///mean temperature,metallicty,star formation rate
    Double_t Temp_mean_gas_sf, Z_mean_gas_sf, SFR_mean_gas_sf;
    //@}

    ///\name gas star forming radial profiles
    //@{
    vector<unsigned int> aperture_npart_gas_sf;
    vector<float> aperture_mass_gas_sf;
    vector<float> aperture_veldisp_gas_sf;
    vector<float> aperture_vrdisp_gas_sf;
    vector<float> aperture_rhalfmass_gas_sf;
    vector<float> aperture_Z_gas_sf;
    vector<Coordinate> aperture_L_gas_sf;
    vector<Coordinate> aperture_mass_proj_gas_sf;
    vector<Coordinate> aperture_rhalfmass_proj_gas_sf;
    vector<Coordinate> aperture_Z_proj_gas_sf;
    vector<unsigned int> profile_npart_gas_sf;
    vector<unsigned int> profile_npart_inclusive_gas_sf;
    vector<float> profile_mass_gas_sf;
    vector<float> profile_mass_inclusive_gas_sf;
    vector<Coordinate> profile_L_gas_sf;
    //@}

    vector<Double_t> SO_mass_gas_sf;
    vector<Coordinate> SO_angularmomentum_gas_sf;

    ///\name star forming gas specific quantities
    //@{
    ///number of particles
    int n_gas_nsf;
    ///mass
    Double_t M_gas_nsf, M_gas_nsf_rvmax,M_gas_nsf_30kpc,M_gas_nsf_50kpc, M_gas_nsf_500c;
    ///mass in spherical overdensities
    Double_t M_200crit_gas_nsf, M_200mean_gas_nsf, M_BN98_gas_nsf;
    ///mass in spherical overdensities inclusive of all masses
    Double_t M_200crit_excl_gas_nsf, M_200mean_excl_gas_nsf, M_BN98_excl_gas_nsf;
    ///velocity/angular momentum info
    Double_t Krot_gas_nsf;
    Coordinate L_gas_nsf;
    ///physical properties for angular momentum (can be inclusive or exclusive )
    Coordinate L_200crit_gas_nsf, L_200mean_gas_nsf, L_BN98_gas_nsf;
    ///physical properties for angular momentum exclusiveto object
    Coordinate L_200crit_excl_gas_nsf, L_200mean_excl_gas_nsf, L_BN98_excl_gas_nsf;
    //dispersion
    Double_t sigV_gas_nsf;
    ///morphology
    Double_t MassTwiceRhalfmass_gas_nsf, Rhalfmass_gas_nsf, q_gas_nsf, s_gas_nsf;
    ///mass weighted sum of temperature, metallicty, star formation rate
    Double_t Temp_gas_nsf, Z_gas_nsf;
    ///mean temperature,metallicty,star formation rate
    Double_t Temp_mean_gas_nsf, Z_mean_gas_nsf;
    //@}

    ///\name gas star forming radial profiles
    //@{
    vector<unsigned int> aperture_npart_gas_nsf;
    vector<float> aperture_mass_gas_nsf;
    vector<float> aperture_veldisp_gas_nsf;
    vector<float> aperture_vrdisp_gas_nsf;
    vector<float> aperture_rhalfmass_gas_nsf;
    vector<float> aperture_Z_gas_nsf;
    vector<Coordinate> aperture_L_gas_nsf;
    vector<Coordinate> aperture_mass_proj_gas_nsf;
    vector<Coordinate> aperture_rhalfmass_proj_gas_nsf;
    vector<Coordinate> aperture_Z_proj_gas_nsf;
    vector<unsigned int> profile_npart_gas_nsf;
    vector<unsigned int> profile_npart_inclusive_gas_nsf;
    vector<float> profile_mass_gas_nsf;
    vector<float> profile_mass_inclusive_gas_nsf;
    vector<Coordinate> profile_L_gas_nsf;
    //@}

    vector<Double_t> SO_mass_gas_nsf;
    vector<Coordinate> SO_angularmomentum_gas_nsf;
#endif
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
    Double_t MassTwiceRhalfmass_star, Rhalfmass_star,q_star,s_star;
    Matrix eigvec_star;
    ///mean age,metallicty
    Double_t t_star,Z_star;
    ///mean age,metallicty
    Double_t t_mean_star,Z_mean_star;
    ///physical properties for dynamical state
    Double_t Efrac_star,Pot_star,T_star;
    //@}

    ///\name stellar radial profiles
    //@{
    vector<unsigned int> aperture_npart_star;
    vector<float> aperture_mass_star;
    vector<float> aperture_veldisp_star;
    vector<float> aperture_vrdisp_star;
    vector<float> aperture_rhalfmass_star;
    vector<float> aperture_Z_star;
    vector<Coordinate> aperture_L_star;
    vector<Coordinate> aperture_mass_proj_star;
    vector<Coordinate> aperture_rhalfmass_proj_star;
    vector<Coordinate> aperture_Z_proj_star;
    vector<unsigned int> profile_npart_star;
    vector<unsigned int> profile_npart_inclusive_star;
    vector<float> profile_mass_star;
    vector<float> profile_mass_inclusive_star;
    vector<Coordinate> profile_L_star;
    //@}

    vector<Double_t> SO_mass_star;
    vector<Coordinate> SO_angularmomentum_star;
#endif

#ifdef BHON
    ///\name black hole specific quantities
    //@{
    ///number of BH
    int n_bh;
    ///mass
    Double_t M_bh, M_bh_mostmassive;
    ///mean accretion rate, metallicty
    Double_t acc_bh, acc_bh_mostmassive;

    ///\name blackhole aperture/radial profiles
    //@{
    vector<int> aperture_npart_bh;
    vector<float> aperture_mass_bh;
    vector<Coordinate> aperture_mass_proj_bh;
    vector<Coordinate> aperture_L_bh;
    //@}

    vector<Double_t> SO_mass_bh;
    vector<Coordinate> SO_angularmomentum_bh;
    //@}
#endif

#ifdef HIGHRES
    ///\name low resolution interloper particle specific quantities
    //@{
    ///number of interloper low res particles
    int n_interloper;
    ///mass
    Double_t M_interloper;
    ///mass in spherical overdensities
    Double_t M_200crit_interloper, M_200mean_interloper, M_BN98_interloper;
    ///mass in spherical overdensities inclusive of all masses
    Double_t M_200crit_excl_interloper, M_200mean_excl_interloper, M_BN98_excl_interloper;

    vector<unsigned int> aperture_npart_interloper;
    vector<float> aperture_mass_interloper;
    vector<Coordinate> aperture_mass_proj_interloper;
    vector<unsigned int> profile_npart_interloper;
    vector<unsigned int> profile_npart_inclusive_interloper;
    vector<float> profile_mass_interloper;
    vector<float> profile_mass_inclusive_interloper;

    vector<Double_t> SO_mass_interloper;
    //@}
#endif

    /// \name extra hydro/star/bh properties such as chemistry/feedback/metal production

    //@{
#if defined(GASON)
    HydroProperties hydroprop;
    vector<HydroProperties> aperture_properties_gas;
#if defined(STARON)
    vector<HydroProperties> aperture_properties_gas_sf;
    vector<HydroProperties> aperture_properties_gas_nsf;
#endif
#endif
#if defined(STARON)
    StarProperties starprop;
    vector<StarProperties> aperture_properties_star;
#endif
#if defined(BHON)
    BHProperties bhprop;
    vector<BHProperties> aperture_properties_bh;
#endif
#if defined(EXTRADMON)
    Int_t n_dm;
    ExtraDMProperties extradmprop;
    vector<ExtraDMProperties> aperture_properties_extra_dm;
#endif
    //@}

    PropData()
    {
        num=gNFOF=gN6DFOF=0;
        gmass=gsize=gRmbp=gmaxvel=gRmaxvel=gRvir=gR200m=gR200c=gRhalfmass=gMassTwiceRhalfmass=Efrac=Pot=T=0.;
        gMFOF=gM6DFOF=0;
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
        MassTwiceRhalfmass_gas=Rhalfmass_gas=0;
        eigvec_gas=Matrix(1,0,0,0,1,0,0,0,1);
        Temp_gas=Z_gas=SFR_gas=0.0;
        Temp_mean_gas=Z_mean_gas=SFR_mean_gas=0.0;
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
#ifdef STARON
        n_gas_sf = n_gas_nsf = 0;
        M_gas_sf=M_gas_sf_rvmax=M_gas_sf_30kpc=M_gas_sf_50kpc=0;
        L_gas_sf[0]=L_gas_sf[1]=L_gas_sf[2]=0;
        q_gas_sf=s_gas_sf=1.0;
        MassTwiceRhalfmass_gas_sf=Rhalfmass_gas_sf=0;
        Temp_gas_sf=Z_gas_sf=0.0;
        Temp_mean_gas_sf=Z_mean_gas_sf=0.0;
        sigV_gas_sf=0;

        M_200mean_gas_sf=M_200crit_gas_sf=M_BN98_gas_sf=0;
        M_200mean_excl_gas_sf=M_200crit_excl_gas_sf=M_BN98_excl_gas_sf=0;
        L_200crit_gas_sf[0]=L_200crit_gas_sf[1]=L_200crit_gas_sf[2]=0;
        L_200mean_gas_sf[0]=L_200mean_gas_sf[1]=L_200mean_gas_sf[2]=0;
        L_BN98_gas_sf[0]=L_BN98_gas_sf[1]=L_BN98_gas_sf[2]=0;
        L_200crit_excl_gas_sf[0]=L_200crit_excl_gas_sf[1]=L_200crit_excl_gas_sf[2]=0;
        L_200mean_excl_gas_sf[0]=L_200mean_excl_gas_sf[1]=L_200mean_excl_gas_sf[2]=0;
        L_BN98_excl_gas_sf[0]=L_BN98_excl_gas_sf[1]=L_BN98_excl_gas_sf[2]=0;

        M_gas_nsf=M_gas_nsf_rvmax=M_gas_nsf_30kpc=M_gas_nsf_50kpc=0;
        L_gas_nsf[0]=L_gas_nsf[1]=L_gas_nsf[2]=0;
        q_gas_nsf=s_gas_nsf=1.0;
        MassTwiceRhalfmass_gas_nsf=Rhalfmass_gas_nsf=0;
        Temp_gas_nsf=Z_gas_nsf=0.0;
        Temp_mean_gas_nsf=Z_mean_gas_nsf=0.0;
        sigV_gas_nsf=0;
        M_200mean_gas_nsf=M_200crit_gas_nsf=M_BN98_gas_nsf=0;
        M_200mean_excl_gas_nsf=M_200crit_excl_gas_nsf=M_BN98_excl_gas_nsf=0;
        L_200crit_gas_nsf[0]=L_200crit_gas_nsf[1]=L_200crit_gas_nsf[2]=0;
        L_200mean_gas_nsf[0]=L_200mean_gas_nsf[1]=L_200mean_gas_nsf[2]=0;
        L_BN98_gas_nsf[0]=L_BN98_gas_nsf[1]=L_BN98_gas_nsf[2]=0;
        L_200crit_excl_gas_nsf[0]=L_200crit_excl_gas_nsf[1]=L_200crit_excl_gas_nsf[2]=0;
        L_200mean_excl_gas_nsf[0]=L_200mean_excl_gas_nsf[1]=L_200mean_excl_gas_nsf[2]=0;
        L_BN98_excl_gas_nsf[0]=L_BN98_excl_gas_nsf[1]=L_BN98_excl_gas_nsf[2]=0;
#endif
#endif
#ifdef STARON
        M_star_rvmax=M_star_30kpc=M_star_50kpc=0;
        n_star=M_star=Efrac_star=0;
        cm_star[0]=cm_star[1]=cm_star[2]=cmvel_star[0]=cmvel_star[1]=cmvel_star[2]=0.;
        L_star[0]=L_star[1]=L_star[2]=0;
        q_star=s_star=1.0;
        MassTwiceRhalfmass_star=Rhalfmass_star=0;
        eigvec_star=Matrix(1,0,0,0,1,0,0,0,1);
        t_star=Z_star=0.;
        t_mean_star=Z_mean_star=0.;
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
        M_bh_mostmassive=0;
        acc_bh=0;
        acc_bh_mostmassive=0;
#endif
#ifdef HIGHRES
        n_interloper=M_interloper=0;
#endif
#ifdef EXTRADMON
        n_dm = 0;
#endif
    }
    ///equals operator, useful if want inclusive information before substructure search
    PropData& operator=(const PropData &p) = default;
    /*
    PropData& operator=(const PropData &p) {
        num=p.num;
        gcm=p.gcm;gcmvel=p.gcmvel;
        gposmbp=p.gposmbp;gvelmbp=p.gvelmbp;
        gposminpot=p.gposminpot;gvelminpot=p.gvelminpot;
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
#ifdef STARON
        M_200mean_gas_sf=p.M_200mean_gas_sf;
        M_200crit_gas_sf=p.M_200crit_gas_sf;
        M_BN98_gas_sf=p.M_BN98_gas_sf;
        M_200mean_excl_gas_sf=p.M_200mean_excl_gas_sf;
        M_200crit_excl_gas_sf=p.M_200crit_excl_gas_sf;
        M_BN98_excl_gas_sf=p.M_BN98_excl_gas_sf;
        L_200mean_gas_sf=p.L_200mean_gas_sf;
        L_200crit_gas_sf=p.L_200crit_gas_sf;
        L_BN98_gas_sf=p.L_BN98_gas_sf;
        L_200mean_excl_gas_sf=p.L_200mean_excl_gas_sf;
        L_200crit_excl_gas_sf=p.L_200crit_excl_gas_sf;
        L_BN98_excl_gas_sf=p.L_BN98_excl_gas_sf;

        M_200mean_gas_nsf=p.M_200mean_gas_nsf;
        M_200crit_gas_nsf=p.M_200crit_gas_nsf;
        M_BN98_gas_nsf=p.M_BN98_gas_nsf;
        M_200mean_excl_gas_nsf=p.M_200mean_excl_gas_nsf;
        M_200crit_excl_gas_nsf=p.M_200crit_excl_gas_nsf;
        M_BN98_excl_gas_nsf=p.M_BN98_excl_gas_nsf;
        L_200mean_gas_nsf=p.L_200mean_gas_nsf;
        L_200crit_gas_nsf=p.L_200crit_gas_nsf;
        L_BN98_gas_nsf=p.L_BN98_gas_nsf;
        L_200mean_excl_gas_nsf=p.L_200mean_excl_gas_nsf;
        L_200crit_excl_gas_nsf=p.L_200crit_excl_gas_nsf;
        L_BN98_excl_gas_nsf=p.L_BN98_excl_gas_nsf;
#endif
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
        aperture_npart=p.aperture_npart;
        aperture_mass=p.aperture_mass;
        aperture_veldisp=p.aperture_veldisp;
        aperture_vrdisp=p.aperture_vrdisp;
        aperture_rhalfmass=p.aperture_rhalfmass;
#if defined(GASON) || defined(STARON) || defined(BHON)
        aperture_npart_dm=p.aperture_npart_dm;
        aperture_mass_dm=p.aperture_mass_dm;
        aperture_veldisp_dm=p.aperture_veldisp_dm;
        aperture_vrdisp_dm=p.aperture_vrdisp_dm;
        aperture_rhalfmass_dm=p.aperture_rhalfmass_dm;
        #endif
#ifdef GASON
        aperture_npart_gas=p.aperture_npart_gas;
        aperture_mass_gas=p.aperture_mass_gas;
        aperture_veldisp_gas=p.aperture_veldisp_gas;
        aperture_rhalfmass_gas=p.aperture_rhalfmass_gas;
#ifdef STARON
        aperture_SFR_gas=p.aperture_SFR_gas;
        aperture_Z_gas=p.aperture_Z_gas;
        aperture_npart_gas_sf=p.aperture_npart_gas_sf;
        aperture_npart_gas_nsf=p.aperture_npart_gas_nsf;
        aperture_mass_gas_sf=p.aperture_mass_gas_sf;
        aperture_mass_gas_nsf=p.aperture_mass_gas_nsf;
        aperture_veldisp_gas_sf=p.aperture_veldisp_gas_sf;
        aperture_veldisp_gas_nsf=p.aperture_veldisp_gas_nsf;
        aperture_vrdisp_gas_sf=p.aperture_vrdisp_gas_sf;
        aperture_vrdisp_gas_nsf=p.aperture_vrdisp_gas_nsf;
        aperture_rhalfmass_gas_sf=p.aperture_rhalfmass_gas_sf;
        aperture_rhalfmass_gas_nsf=p.aperture_rhalfmass_gas_nsf;
        aperture_Z_gas_sf=p.aperture_Z_gas_sf;
        aperture_Z_gas_nsf=p.aperture_Z_gas_nsf;
#endif
#endif
#ifdef STARON
        aperture_npart_star=p.aperture_npart_star;
        aperture_mass_star=p.aperture_mass_star;
        aperture_veldisp_star=p.aperture_veldisp_star;
        aperture_vrdisp_star=p.aperture_vrdisp_star;
        aperture_rhalfmass_star=p.aperture_rhalfmass_star;
        aperture_Z_star=p.aperture_Z_star;
#endif
        aperture_mass_proj=p.aperture_mass_proj;
        aperture_rhalfmass_proj=p.aperture_rhalfmass_proj;
#ifdef GASON
        aperture_mass_proj_gas=p.aperture_mass_proj_gas;
        aperture_rhalfmass_proj_gas=p.aperture_rhalfmass_proj_gas;
#ifdef STARON
        aperture_SFR_proj_gas=p.aperture_SFR_proj_gas;
        aperture_Z_proj_gas=p.aperture_Z_proj_gas;
        aperture_mass_proj_gas_sf=p.aperture_mass_proj_gas_sf;
        aperture_mass_proj_gas_nsf=p.aperture_mass_proj_gas_nsf;
        aperture_rhalfmass_proj_gas_sf=p.aperture_rhalfmass_proj_gas_sf;
        aperture_rhalfmass_proj_gas_nsf=p.aperture_rhalfmass_proj_gas_nsf;
        aperture_Z_proj_gas_sf=p.aperture_Z_proj_gas_sf;
        aperture_Z_proj_gas_nsf=p.aperture_Z_proj_gas_nsf;
#endif
#endif
#ifdef STARON
        aperture_mass_proj_star=p.aperture_mass_proj_star;
        aperture_rhalfmass_proj_star=p.aperture_rhalfmass_proj_star;
        aperture_Z_proj_star=p.aperture_Z_proj_star;
#endif
        profile_npart=p.profile_npart;
        profile_mass=p.profile_mass;
        profile_npart_inclusive=p.profile_npart_inclusive;
        profile_mass_inclusive=p.profile_mass_inclusive;
#ifdef GASON
        profile_npart_gas=p.profile_npart_gas;
        profile_mass_gas=p.profile_mass_gas;
        profile_npart_inclusive_gas=p.profile_npart_inclusive_gas;
        profile_mass_inclusive_gas=p.profile_mass_inclusive_gas;
#ifdef STARON
        profile_npart_gas_sf=p.profile_npart_gas_sf;
        profile_mass_gas_sf=p.profile_mass_gas_sf;
        profile_npart_inclusive_gas_sf=p.profile_npart_inclusive_gas_sf;
        profile_mass_inclusive_gas_sf=p.profile_mass_inclusive_gas_sf;
        profile_npart_gas_nsf=p.profile_npart_gas_nsf;
        profile_mass_gas_nsf=p.profile_mass_gas_nsf;
        profile_npart_inclusive_gas_nsf=p.profile_npart_inclusive_gas_nsf;
        profile_mass_inclusive_gas_nsf=p.profile_mass_inclusive_gas_nsf;
#endif
#endif
#ifdef STARON
        profile_npart_star=p.profile_npart_star;
        profile_mass_star=p.profile_mass_star;
        profile_npart_inclusive_star=p.profile_npart_inclusive_star;
        profile_mass_inclusive_star=p.profile_mass_inclusive_star;
#endif
        return *this;
    }
    */

    //allocate memory for profiles
    void Allocate(Options &opt) {
        AllocateApertures(opt);
        AllocateProfiles(opt);
        AllocateSOs(opt);
    }
    void AllocateApertures(Options &opt)
    {
        if (opt.iaperturecalc && opt.aperturenum>0) {
            aperture_npart.resize(opt.aperturenum);
            aperture_mass.resize(opt.aperturenum);
            aperture_veldisp.resize(opt.aperturenum);
            aperture_vrdisp.resize(opt.aperturenum);
            aperture_rhalfmass.resize(opt.aperturenum);
#ifdef GASON
            aperture_npart_gas.resize(opt.aperturenum);
            aperture_mass_gas.resize(opt.aperturenum);
            aperture_veldisp_gas.resize(opt.aperturenum);
            aperture_vrdisp_gas.resize(opt.aperturenum);
            aperture_rhalfmass_gas.resize(opt.aperturenum);
            if (opt.gas_extraprop_aperture_calc) aperture_properties_gas.resize(opt.aperturenum);
#ifdef STARON
            aperture_SFR_gas.resize(opt.aperturenum);
            aperture_Z_gas.resize(opt.aperturenum);
            aperture_npart_gas_sf.resize(opt.aperturenum);
            aperture_npart_gas_nsf.resize(opt.aperturenum);
            aperture_mass_gas_sf.resize(opt.aperturenum);
            aperture_mass_gas_nsf.resize(opt.aperturenum);
            aperture_veldisp_gas_sf.resize(opt.aperturenum);
            aperture_veldisp_gas_nsf.resize(opt.aperturenum);
            aperture_vrdisp_gas_sf.resize(opt.aperturenum);
            aperture_vrdisp_gas_nsf.resize(opt.aperturenum);
            aperture_rhalfmass_gas_sf.resize(opt.aperturenum);
            aperture_rhalfmass_gas_nsf.resize(opt.aperturenum);
            aperture_Z_gas_sf.resize(opt.aperturenum);
            aperture_Z_gas_nsf.resize(opt.aperturenum);
            // if (opt.gas_extraprop_aperture_calc) aperture_properties_gas_sf.resize(opt.aperturenum);
            // if (opt.gas_extraprop_aperture_calc) aperture_properties_gas_nsf.resize(opt.aperturenum);
#endif
#endif
#ifdef STARON
            aperture_npart_star.resize(opt.aperturenum);
            aperture_mass_star.resize(opt.aperturenum);
            aperture_veldisp_star.resize(opt.aperturenum);
            aperture_vrdisp_star.resize(opt.aperturenum);
            aperture_rhalfmass_star.resize(opt.aperturenum);
            aperture_Z_star.resize(opt.aperturenum);
            if (opt.star_extraprop_aperture_calc) aperture_properties_star.resize(opt.aperturenum);
#endif
#ifdef BHON
            aperture_npart_bh.resize(opt.aperturenum);
            aperture_mass_bh.resize(opt.aperturenum);
            if (opt.bh_extraprop_aperture_calc) aperture_properties_bh.resize(opt.aperturenum);
#endif
#ifdef HIGHRES
            aperture_npart_interloper.resize(opt.aperturenum);
            aperture_mass_interloper.resize(opt.aperturenum);
#endif
#ifdef EXTRADMON
            if (opt.extra_dm_extraprop_aperture_calc) aperture_properties_extra_dm.resize(opt.aperturenum);
#endif
#if defined(GASON) || defined(STARON) || defined(BHON)
            //if searching all types, also store dm only aperture quantities
            if (opt.partsearchtype==PSTALL) {
                aperture_npart_dm.resize(opt.aperturenum);
                aperture_mass_dm.resize(opt.aperturenum);
                aperture_veldisp_dm.resize(opt.aperturenum);
                aperture_vrdisp_dm.resize(opt.aperturenum);
                aperture_rhalfmass_dm.resize(opt.aperturenum);
            }
#endif
            for (auto &x:aperture_npart) x=0;
            for (auto &x:aperture_mass) x=-1;
            for (auto &x:aperture_veldisp) x=0;
            for (auto &x:aperture_rhalfmass) x=-1;
#ifdef GASON
            for (auto &x:aperture_npart_gas) x=0;
            for (auto &x:aperture_mass_gas) x=-1;
            for (auto &x:aperture_veldisp_gas) x=0;
            for (auto &x:aperture_rhalfmass_gas) x=-1;
#ifdef STARON
            for (auto &x:aperture_SFR_gas) x=0;
            for (auto &x:aperture_Z_gas) x=0;
            for (auto &x:aperture_npart_gas_sf) x=0;
            for (auto &x:aperture_mass_gas_sf) x=-1;
            for (auto &x:aperture_npart_gas_nsf) x=0;
            for (auto &x:aperture_mass_gas_nsf) x=-1;
            for (auto &x:aperture_veldisp_gas_sf) x=0;
            for (auto &x:aperture_veldisp_gas_nsf) x=0;
            for (auto &x:aperture_rhalfmass_gas_sf) x=-1;
            for (auto &x:aperture_rhalfmass_gas_nsf) x=-1;
            for (auto &x:aperture_Z_gas_sf) x=0;
            for (auto &x:aperture_Z_gas_nsf) x=0;
#endif
#endif
#ifdef STARON
            for (auto &x:aperture_npart_star) x=0;
            for (auto &x:aperture_mass_star) x=-1;
            for (auto &x:aperture_veldisp_star) x=0;
            for (auto &x:aperture_rhalfmass_star) x=-1;
            for (auto &x:aperture_Z_star) x=0;
#endif
#ifdef HIGHRES
            for (auto &x:aperture_npart_interloper) x=0;
            for (auto &x:aperture_mass_interloper) x=-1;
#endif
#if defined(GASON) || defined(STARON) || defined(BHON)
            if (opt.partsearchtype==PSTALL) {
                for (auto &x:aperture_npart_dm) x=0;
                for (auto &x:aperture_mass_dm) x=-1;
                for (auto &x:aperture_veldisp_dm) x=0;
                for (auto &x:aperture_rhalfmass_dm) x=0;
            }
#endif
        }

        if (opt.iaperturecalc && opt.apertureprojnum>0) {
            aperture_mass_proj.resize(opt.apertureprojnum);
            aperture_rhalfmass_proj.resize(opt.apertureprojnum);
#ifdef GASON
            aperture_mass_proj_gas.resize(opt.apertureprojnum);
            aperture_rhalfmass_proj_gas.resize(opt.apertureprojnum);
#ifdef STARON
            aperture_SFR_proj_gas.resize(opt.apertureprojnum);
            aperture_Z_proj_gas.resize(opt.apertureprojnum);
            aperture_mass_proj_gas_sf.resize(opt.apertureprojnum);
            aperture_mass_proj_gas_nsf.resize(opt.apertureprojnum);
            aperture_rhalfmass_proj_gas_sf.resize(opt.apertureprojnum);
            aperture_rhalfmass_proj_gas_nsf.resize(opt.apertureprojnum);
            aperture_Z_proj_gas_sf.resize(opt.apertureprojnum);
            aperture_Z_proj_gas_nsf.resize(opt.apertureprojnum);
#endif
#endif
#ifdef STARON
            aperture_mass_proj_star.resize(opt.apertureprojnum);
            aperture_rhalfmass_proj_star.resize(opt.apertureprojnum);
            aperture_Z_proj_star.resize(opt.apertureprojnum);
#endif
#ifdef BHON
            aperture_mass_proj_bh.resize(opt.apertureprojnum);
#endif
#ifdef HIGHRES
            aperture_mass_proj_interloper.resize(opt.apertureprojnum);
#endif
            for (auto &x:aperture_mass_proj) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_rhalfmass_proj) x[0]=x[1]=x[2]=-1;
#ifdef GASON
            for (auto &x:aperture_mass_proj_gas) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_rhalfmass_proj_gas) x[0]=x[1]=x[2]=-1;
#ifdef STARON
            for (auto &x:aperture_SFR_proj_gas) x[0]=x[1]=x[2]=0;
            for (auto &x:aperture_Z_proj_gas) x[0]=x[1]=x[2]=0;
            for (auto &x:aperture_mass_proj_gas_sf) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_rhalfmass_proj_gas_sf) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_mass_proj_gas_nsf) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_rhalfmass_proj_gas_nsf) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_Z_proj_gas_sf) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_Z_proj_gas_nsf) x[0]=x[1]=x[2]=-1;
#endif
#endif
#ifdef STARON
            for (auto &x:aperture_mass_proj_star) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_rhalfmass_proj_star) x[0]=x[1]=x[2]=-1;
            for (auto &x:aperture_Z_proj_star) x[0]=x[1]=x[2]=-1;
#endif
#ifdef BHON
            for (auto &x:aperture_mass_proj_bh) x[0]=x[1]=x[2]=-1;
#endif
#ifdef HIGHRES
            for (auto &x:aperture_mass_proj_interloper) x[0]=x[1]=x[2]=-1;
#endif
        }
    }
    void AllocateProfiles(Options &opt)
    {
        if (opt.iprofilecalc && gNFOF>=opt.profileminFOFsize && num>=opt.profileminsize) {
            profile_npart.resize(opt.profilenbins);
            profile_mass.resize(opt.profilenbins);
            for (auto i=0;i<opt.profilenbins;i++) profile_npart[i]=profile_mass[i]=0;
#ifdef GASON
            profile_npart_gas.resize(opt.profilenbins);
            profile_mass_gas.resize(opt.profilenbins);
#ifdef STARON
            profile_npart_gas_sf.resize(opt.profilenbins);
            profile_mass_gas_sf.resize(opt.profilenbins);
            profile_npart_gas_nsf.resize(opt.profilenbins);
            profile_mass_gas_nsf.resize(opt.profilenbins);
#endif
            for (auto i=0;i<opt.profilenbins;i++) profile_npart_gas[i]=profile_mass_gas[i]=0;
#ifdef STARON
            for (auto i=0;i<opt.profilenbins;i++) profile_npart_gas_sf[i]=profile_mass_gas_sf[i]=profile_npart_gas_nsf[i]=profile_mass_gas_nsf[i]=0;
#endif
#endif
#ifdef STARON
            profile_npart_star.resize(opt.profilenbins);
            profile_mass_star.resize(opt.profilenbins);
            for (auto i=0;i<opt.profilenbins;i++) profile_npart_star[i]=profile_mass_star[i]=0;
#endif
            if (opt.iInclusiveHalo>0) {
                profile_npart_inclusive.resize(opt.profilenbins);
                profile_mass_inclusive.resize(opt.profilenbins);
                for (auto i=0;i<opt.profilenbins;i++) profile_npart_inclusive[i]=profile_mass_inclusive[i]=0;
#ifdef GASON
                profile_npart_inclusive_gas.resize(opt.profilenbins);
                profile_mass_inclusive_gas.resize(opt.profilenbins);
#ifdef STARON
                profile_npart_inclusive_gas_sf.resize(opt.profilenbins);
                profile_mass_inclusive_gas_sf.resize(opt.profilenbins);
                profile_npart_inclusive_gas_nsf.resize(opt.profilenbins);
                profile_mass_inclusive_gas_nsf.resize(opt.profilenbins);
#endif
                for (auto i=0;i<opt.profilenbins;i++) profile_npart_inclusive_gas[i]=profile_mass_inclusive_gas[i]=0;
#ifdef STARON
                for (auto i=0;i<opt.profilenbins;i++) profile_npart_inclusive_gas_sf[i]=profile_mass_inclusive_gas_sf[i]=profile_npart_inclusive_gas_nsf[i]=profile_mass_inclusive_gas_nsf[i]=0;
#endif
#endif
#ifdef STARON
                profile_npart_inclusive_star.resize(opt.profilenbins);
                profile_mass_inclusive_star.resize(opt.profilenbins);
                for (auto i=0;i<opt.profilenbins;i++) profile_npart_inclusive_star[i]=profile_mass_inclusive_star[i]=0;
#endif
            }

        }
    }
    void AllocateSOs(Options &opt)
    {
        if (opt.SOnum>0) {
            SO_mass.resize(opt.SOnum);
            SO_radius.resize(opt.SOnum);
            for (auto &x:SO_mass) x=0;
            for (auto &x:SO_radius) x=0;
            if (opt.iextrahalooutput) {
                SO_angularmomentum.resize(opt.SOnum);
                for (auto &x:SO_angularmomentum) {x[0]=x[1]=x[2]=0;}
#ifdef GASON
                if (opt.iextragasoutput) {
                    SO_mass_gas.resize(opt.SOnum);
                    for (auto &x:SO_mass_gas) x=0;
                    SO_angularmomentum_gas.resize(opt.SOnum);
                    for (auto &x:SO_angularmomentum_gas) {x[0]=x[1]=x[2]=0;}
#ifdef STARON
#endif
                }
#endif
#ifdef STARON
                if (opt.iextrastaroutput) {
                    SO_mass_star.resize(opt.SOnum);
                    for (auto &x:SO_mass_star) x=0;
                    SO_angularmomentum_star.resize(opt.SOnum);
                    for (auto &x:SO_angularmomentum_star) {x[0]=x[1]=x[2]=0;}
                }
#endif
#ifdef HIGHRES
                if (opt.iextrainterloperoutput) {
                    SO_mass_interloper.resize(opt.SOnum);
                    for (auto &x:SO_mass_interloper) x=0;
                }
#endif
            }
        }
    }
    void CopyProfileToInclusive(Options &opt) {
        for (auto i=0;i<opt.profilenbins;i++) {
            profile_npart_inclusive[i]=profile_npart[i];
            profile_mass_inclusive[i]=profile_mass[i];
            profile_npart[i]=profile_mass[i]=0;
#ifdef GASON
            profile_npart_inclusive_gas[i]=profile_npart_gas[i];
            profile_mass_inclusive_gas[i]=profile_mass_gas[i];
            profile_npart_gas[i]=profile_mass_gas[i]=0;
#ifdef STARON
            profile_npart_inclusive_gas_sf[i]=profile_npart_gas_sf[i];
            profile_mass_inclusive_gas_sf[i]=profile_mass_gas_sf[i];
            profile_npart_gas_sf[i]=profile_mass_gas_sf[i]=0;
            profile_npart_inclusive_gas_nsf[i]=profile_npart_gas_nsf[i];
            profile_mass_inclusive_gas_nsf[i]=profile_mass_gas_nsf[i];
            profile_npart_gas_nsf[i]=profile_mass_gas_nsf[i]=0;
#endif
#endif
#ifdef STARON
            profile_npart_inclusive_star[i]=profile_npart_star[i];
            profile_mass_inclusive_star[i]=profile_mass_star[i];
            profile_npart_star[i]=profile_mass_star[i]=0;
#endif
        }
    }

    ///converts the properties data into comoving little h values
    ///so masses, positions have little h values and positions are comoving
    void ConverttoComove(Options &opt){
        gcm=gcm*opt.h/opt.a;
        gposmbp=gposmbp*opt.h/opt.a;
        gposminpot=gposminpot*opt.h/opt.a;
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
        gMassTwiceRhalfmass*=opt.h;
        gRhalfmass*=opt.h/opt.a;
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
        L_gas=L_gas*opt.h*opt.h/opt.a;
        MassTwiceRhalfmass_gas*=opt.h;
        Rhalfmass_gas*=opt.h/opt.a;

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
        #ifdef STARON
        M_gas_sf*=opt.h;
        L_gas_sf*=opt.h*opt.h/opt.a;
        MassTwiceRhalfmass_gas_sf*=opt.h;
        Rhalfmass_gas_sf*=opt.h/opt.a;

        M_gas_nsf*=opt.h;
        L_gas_nsf*=opt.h*opt.h/opt.a;
        MassTwiceRhalfmass_gas_nsf*=opt.h;
        Rhalfmass_gas_nsf*=opt.h/opt.a;

        if (opt.iextragasoutput) {
            M_200mean_gas_sf*=opt.h;
            M_200crit_gas_sf*=opt.h;
            M_BN98_gas_sf*=opt.h;
            M_200mean_excl_gas_sf*=opt.h;
            M_200crit_excl_gas_sf*=opt.h;
            M_BN98_excl_gas_sf*=opt.h;
            L_200crit_gas_sf*=opt.h*opt.h/opt.a;
            L_200mean_gas_sf*=opt.h*opt.h/opt.a;
            L_BN98_gas_sf*=opt.h*opt.h/opt.a;
            L_200crit_excl_gas_sf*=opt.h*opt.h/opt.a;
            L_200mean_excl_gas_sf*=opt.h*opt.h/opt.a;
            L_BN98_excl_gas_sf*=opt.h*opt.h/opt.a;

            M_200mean_gas_nsf*=opt.h;
            M_200crit_gas_nsf*=opt.h;
            M_BN98_gas_nsf*=opt.h;
            M_200mean_excl_gas_nsf*=opt.h;
            M_200crit_excl_gas_nsf*=opt.h;
            M_BN98_excl_gas_nsf*=opt.h;
            L_200crit_gas_nsf*=opt.h*opt.h/opt.a;
            L_200mean_gas_nsf*=opt.h*opt.h/opt.a;
            L_BN98_gas_nsf*=opt.h*opt.h/opt.a;
            L_200crit_excl_gas_nsf*=opt.h*opt.h/opt.a;
            L_200mean_excl_gas_nsf*=opt.h*opt.h/opt.a;
            L_BN98_excl_gas_nsf*=opt.h*opt.h/opt.a;
        }
        #endif

#endif
#ifdef STARON
        M_star*=opt.h;
        M_star_rvmax*=opt.h;
        M_star_30kpc*=opt.h;
        M_star_50kpc*=opt.h;
        M_star_500c*=opt.h;
        cm_star=cm_star*opt.h/opt.a;
        MassTwiceRhalfmass_star*=opt.h/opt.a;
        Rhalfmass_star*=opt.h/opt.a;
        L_star=L_star*opt.h*opt.h/opt.a;
#endif
#ifdef BHON
        M_bh*=opt.h;
#endif
#ifdef HIGHRES
        M_interloper*=opt.h;
#endif
        if (opt.iaperturecalc) {
            for (auto i=0;i<opt.iaperturecalc;i++) {
                aperture_mass[i]*=opt.h;
#ifdef GASON
                aperture_mass_gas[i]*=opt.h;
#ifdef STARON
                aperture_mass_gas_sf[i]*=opt.h;
                aperture_mass_gas_nsf[i]*=opt.h;
#endif
#endif
#ifdef STARON
                aperture_mass_star[i]*=opt.h;
#endif
            }
        }
        if (opt.SOnum>0) {
            for (auto i=0;i<opt.SOnum;i++) {
                SO_mass[i] *= opt.h;
                SO_radius[i] *= opt.h/opt.a;
                if (opt.iextrahalooutput) {
                    SO_angularmomentum[i]*=(opt.h*opt.h/opt.a);
#ifdef GASON
                    SO_mass_gas[i] *= opt.h;
                    SO_angularmomentum_gas[i]*=(opt.h*opt.h/opt.a);
#ifdef STARON
#endif
#endif
#ifdef STARON
                    SO_mass_star[i] *= opt.h;
                    SO_angularmomentum_star[i]*=(opt.h*opt.h/opt.a);
#endif
                }
            }
        }
    }

    void ConvertProfilestoComove(Options &opt){
        for (auto i=0;i<opt.profilenbins;i++) {
            profile_mass[i]*=opt.h;
#ifdef GASON
            profile_mass_gas[i]*=opt.h;
#ifdef STARON
            profile_mass_gas_sf[i]*=opt.h;
            profile_mass_gas_nsf[i]*=opt.h;
#endif
#endif
#ifdef STARON
            profile_mass_star[i]*=opt.h;
#endif
            }
        if (opt.iInclusiveHalo) {
            for (auto i=0;i<opt.profilenbins;i++) {
                profile_mass_inclusive[i]*=opt.h;
#ifdef GASON
                profile_mass_inclusive_gas[i]*=opt.h;
#ifdef STARON
                profile_mass_inclusive_gas_sf[i]*=opt.h;
                profile_mass_inclusive_gas_nsf[i]*=opt.h;
#endif
#endif
#ifdef STARON
                profile_mass_inclusive_star[i]*=opt.h;
#endif
            }
        }
    }

    void WriteBinary(fstream &Fout, Options &opt);
    void WriteAscii(fstream &Fout, Options &opt);
    ///write (append) the properties data to an already open ascii file
#ifdef USEHDF
    ///write (append) the properties data to an already open hdf file
    //void WriteHDF(H5File &Fhdf, DataSpace *&dataspaces, DataSet *&datasets, Options&opt){};
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
    // vector<PredType> predtypeinfo;
    vector<hid_t> hdfpredtypeinfo;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospredtypeinfo;
#endif
    PropDataHeader(Options&opt);
};


/*! Structures stores profile info of the data writen by the \ref PropData profiles data structures,
    specifically the \ref PropData::WriteProfileBinary, \ref PropData::WriteProfileAscii, \ref PropData::WriteProfileHDF routines
*/

struct ProfileDataHeader{
    //list the header info
    vector<string> headerdatainfo;
#ifdef USEHDF
    vector<hid_t> hdfpredtypeinfo;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospredtypeinfo;
#endif
    int numberscalarentries, numberarrayallgroupentries, numberarrayhaloentries;
    int offsetscalarentries, offsetarrayallgroupentries, offsetarrayhaloentries;

    ProfileDataHeader(Options&opt){
        int sizeval;
#ifdef USEHDF
        vector<hid_t> hdfdesiredproprealtype;
        if (sizeof(Double_t)==sizeof(double)) hdfdesiredproprealtype.push_back(H5T_NATIVE_DOUBLE);
        else hdfdesiredproprealtype.push_back(H5T_NATIVE_FLOAT);
#endif
#ifdef USEADIOS
        vector<ADIOS_DATATYPES> desiredadiosproprealtype;
        if (sizeof(Double_t)==sizeof(double)) desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_double);
        else desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_real);
#endif

        offsetscalarentries=0;
        headerdatainfo.push_back("ID");
#ifdef USEHDF
        hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
        //if normalisation is phys then no need for writing normalisation block
        if (opt.iprofilenorm != PROFILERNORMPHYS) {
        headerdatainfo.push_back(opt.profileradnormstring);
#ifdef USEHDF
        hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
        adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
        numberscalarentries=headerdatainfo.size();

	    offsetarrayallgroupentries=headerdatainfo.size();
        headerdatainfo.push_back("Npart_profile");
#ifdef GASON
        headerdatainfo.push_back("Npart_profile_gas");
#ifdef STARON
        headerdatainfo.push_back("Npart_profile_gas_sf");
        headerdatainfo.push_back("Npart_profile_gas_nsf");
#endif
#endif
#ifdef STARON
        headerdatainfo.push_back("Npart_profile_star");
#endif
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(PredType::STD_U32LE);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(H5T_NATIVE_UINT);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
#endif

        headerdatainfo.push_back("Mass_profile");
#ifdef GASON
        headerdatainfo.push_back("Mass_profile_gas");
#ifdef STARON
        headerdatainfo.push_back("Mass_profile_gas_sf");
        headerdatainfo.push_back("Mass_profile_gas_nsf");
#endif
#endif
#ifdef STARON
        headerdatainfo.push_back("Mass_profile_star");
#endif

#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(PredType::NATIVE_FLOAT);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(H5T_NATIVE_FLOAT);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_real);
#endif
        numberarrayallgroupentries=headerdatainfo.size()-offsetarrayallgroupentries;

        //stuff for inclusive halo/SO profiles
        if (opt.iInclusiveHalo >0) {
        offsetarrayhaloentries=headerdatainfo.size();
        headerdatainfo.push_back("Npart_inclusive_profile");
#ifdef GASON
        headerdatainfo.push_back("Npart_inclusive_profile_gas");
#ifdef STARON
        headerdatainfo.push_back("Npart_inclusive_profile_gas_sf");
        headerdatainfo.push_back("Npart_inclusive_profile_gas_nsf");
#endif
#endif
#ifdef STARON
        headerdatainfo.push_back("Npart_inclusive_profile_star");
#endif
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(PredType::STD_U32LE);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(H5T_NATIVE_UINT);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
#endif

        headerdatainfo.push_back("Mass_inclusive_profile");
#ifdef GASON
        headerdatainfo.push_back("Mass_inclusive_profile_gas");
#ifdef STARON
        headerdatainfo.push_back("Mass_inclusive_profile_gas_sf");
        headerdatainfo.push_back("Mass_inclusive_profile_gas_nsf");
#endif
#endif
#ifdef STARON
        headerdatainfo.push_back("Mass_inclusive_profile_star");
#endif
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(PredType::NATIVE_FLOAT);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(H5T_NATIVE_FLOAT);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_real);
#endif
        numberarrayhaloentries=headerdatainfo.size()-offsetarrayhaloentries;
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
        Pparenthead=new Particle*[numgroups+1];
        gidhead=new Int_t*[numgroups+1];
        gidparenthead=new Int_t*[numgroups+1];
        giduberparenthead=new Int_t*[numgroups+1];
        stypeinlevel=new Int_t[numgroups+1];
        nextlevel=NULL;
    }
    ///initialize
    void Initialize(){
        for (Int_t i=1;i<=nsinlevel;i++) {gidhead[i]=NULL;gidparenthead[i]=NULL;giduberparenthead[i]=NULL;}
    }
    ~StrucLevelData(){
        if (nextlevel!=NULL) delete nextlevel;
        nextlevel=NULL;
        delete[] Phead;
        delete[] Pparenthead;
        delete[] gidhead;
        delete[] gidparenthead;
        delete[] giduberparenthead;
        delete[] stypeinlevel;
    }
};

#if defined(USEHDF)||defined(USEADIOS)
///store the names of datasets in catalog output
struct DataGroupNames {
    ///store names of catalog group files
    vector<string> prop;
#ifdef USEHDF
    //store the data type
    // vector<PredType> propdatatype;
    vector<hid_t> hdfpropdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospropdatatype;
#endif

    ///store names of catalog group files
    vector<string> group;
#ifdef USEHDF
    // vector<PredType> groupdatatype;
    vector<hid_t> hdfgroupdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiosgroupdatatype;
#endif

    ///store the names of catalog particle files
    vector<string> part;
#ifdef USEHDF
    // vector<PredType> partdatatype;
    vector<hid_t> hdfpartdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiospartdatatype;
#endif

    ///store the names of catalog particle type files
    vector<string> types;
#ifdef USEHDF
    // vector<PredType> typesdatatype;
    vector<hid_t> hdftypesdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiostypesdatatype;
#endif

    ///store the names of hierarchy files
    vector<string> hierarchy;
#ifdef USEHDF
    // vector<PredType> hierarchydatatype;
    vector<hid_t> hdfhierarchydatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adioshierarchydatatype;
#endif

    ///store names of SO files
    vector<string> SO;
#ifdef USEHDF
    // vector<PredType> SOdatatype;
    vector<hid_t> hdfSOdatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> SOdatatype;
#endif

    //store names of profile files
    vector<string> profile;
#ifdef USEHDF
    //store the data type
    // vector<PredType> profiledatatype;
    vector<hid_t> hdfprofiledatatype;
#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> adiosprofiledatatype;
#endif

    DataGroupNames(){
#ifdef USEHDF
        vector<hid_t> hdfdesiredproprealtype;
        if (sizeof(Double_t)==sizeof(double)) hdfdesiredproprealtype.push_back(H5T_NATIVE_DOUBLE);
        else hdfdesiredproprealtype.push_back(H5T_NATIVE_FLOAT);
#endif
#ifdef USEADIOS
        vector<ADIOS_DATATYPES> desiredadiosproprealtype;
        if (sizeof(Double_t)==sizeof(double)) desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_double);
        else desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_real);
#endif
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
#if defined(GASON) || defined(STARON) || defined(BHON)
        prop.push_back("Metallicity_unit_to_solar");
        prop.push_back("SFR_unit_to_solarmassperyear");
        prop.push_back("Stellar_age_unit_to_yr");
#endif
#ifdef USEHDF
        hdfpropdatatype.push_back(H5T_NATIVE_INT);
        hdfpropdatatype.push_back(H5T_NATIVE_INT);
        hdfpropdatatype.push_back(H5T_NATIVE_ULONG);
        hdfpropdatatype.push_back(H5T_NATIVE_ULONG);
        hdfpropdatatype.push_back(H5T_NATIVE_UINT);
        hdfpropdatatype.push_back(H5T_NATIVE_UINT);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
#if defined(GASON) || defined(STARON) || defined(BHON)
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
        hdfpropdatatype.push_back(hdfdesiredproprealtype[0]);
#endif

#endif
#ifdef USEADIOS
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adiospropdatatype.push_back(ADIOS_DATATYPES::adios_unsigned_integer);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
#if defined(GASON) || defined(STARON) || defined(BHON)
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
        adiospropdatatype.push_back(desiredadiosproprealtype[0]);
#endif
#endif

        group.push_back("File_id");
        group.push_back("Num_of_files");
        group.push_back("Num_of_groups");
        group.push_back("Total_num_of_groups");
        group.push_back("Group_Size");
        group.push_back("Offset");
        group.push_back("Offset_unbound");
#ifdef USEHDF
        hdfgroupdatatype.push_back(H5T_NATIVE_INT);
        hdfgroupdatatype.push_back(H5T_NATIVE_INT);
        hdfgroupdatatype.push_back(H5T_NATIVE_ULONG);
        hdfgroupdatatype.push_back(H5T_NATIVE_ULONG);
        hdfgroupdatatype.push_back(H5T_NATIVE_ULONG);
        hdfgroupdatatype.push_back(H5T_NATIVE_ULONG);
        hdfgroupdatatype.push_back(H5T_NATIVE_ULONG);

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
        hdfpartdatatype.push_back(H5T_NATIVE_INT);
        hdfpartdatatype.push_back(H5T_NATIVE_INT);
        hdfpartdatatype.push_back(H5T_NATIVE_ULONG);
        hdfpartdatatype.push_back(H5T_NATIVE_ULONG);
        hdfpartdatatype.push_back(H5T_NATIVE_LONG);

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

        hdftypesdatatype.push_back(H5T_NATIVE_INT);
        hdftypesdatatype.push_back(H5T_NATIVE_INT);
        hdftypesdatatype.push_back(H5T_NATIVE_ULONG);
        hdftypesdatatype.push_back(H5T_NATIVE_ULONG);
        hdftypesdatatype.push_back(H5T_NATIVE_USHORT);
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
        hdfhierarchydatatype.push_back(H5T_NATIVE_INT);
        hdfhierarchydatatype.push_back(H5T_NATIVE_INT);
        hdfhierarchydatatype.push_back(H5T_NATIVE_ULONG);
        hdfhierarchydatatype.push_back(H5T_NATIVE_ULONG);
        hdfhierarchydatatype.push_back(H5T_NATIVE_UINT);
        hdfhierarchydatatype.push_back(H5T_NATIVE_LONG);

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
#if defined(GASON) || defined(STARON) || defined(BHON)
        SO.push_back("Particle_types");
#endif

#ifdef USEHDF
        hdfSOdatatype.push_back(H5T_NATIVE_INT);
        hdfSOdatatype.push_back(H5T_NATIVE_INT);
        hdfSOdatatype.push_back(H5T_NATIVE_ULONG);
        hdfSOdatatype.push_back(H5T_NATIVE_ULONG);
        hdfSOdatatype.push_back(H5T_NATIVE_ULONG);
        hdfSOdatatype.push_back(H5T_NATIVE_ULONG);
        hdfSOdatatype.push_back(H5T_NATIVE_UINT);
        hdfSOdatatype.push_back(H5T_NATIVE_ULONG);
        hdfSOdatatype.push_back(H5T_NATIVE_LONG);
        #if defined(GASON) || defined(STARON) || defined(BHON)
        hdfSOdatatype.push_back(H5T_NATIVE_INT);
        #endif

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
#if defined(GASON) || defined(STARON) || defined(BHON)
        adiosSOdatatype.push_back(ADIOS_DATATYPES::adios_integer);
#endif
#endif

        profile.push_back("File_id");
        profile.push_back("Num_of_files");
        profile.push_back("Num_of_groups");
        profile.push_back("Total_num_of_groups");
        profile.push_back("Num_of_halos");
        profile.push_back("Total_num_of_halos");
        profile.push_back("Radial_norm");
        profile.push_back("Inclusive_profiles_flag");
        profile.push_back("Num_of_bin_edges");
        profile.push_back("Radial_bin_edges");
#ifdef USEHDF
        hdfprofiledatatype.push_back(H5T_NATIVE_INT);
        hdfprofiledatatype.push_back(H5T_NATIVE_INT);
        hdfprofiledatatype.push_back(H5T_NATIVE_ULONG);
        hdfprofiledatatype.push_back(H5T_NATIVE_ULONG);
        hdfprofiledatatype.push_back(H5T_NATIVE_ULONG);
        hdfprofiledatatype.push_back(H5T_NATIVE_ULONG);
        hdfprofiledatatype.push_back(H5T_C_S1);
        hdfprofiledatatype.push_back(H5T_NATIVE_INT);
        hdfprofiledatatype.push_back(H5T_NATIVE_INT);
        hdfprofiledatatype.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_unsigned_long);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_string);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosprofiledatatype.push_back(ADIOS_DATATYPES::adios_integer);
        adiosprofiledatatype.push_back(desiredadiosproprealtype[0]);
#endif

    }
};
#endif

///Useful structore to store information of leaf nodes in the tree
struct leaf_node_info{
    int num, numtot;
    Int_t id, istart, iend;
    Coordinate cm;
    Double_t size;
#ifdef USEMPI
    Double_t searchdist;
#endif
};

///if using MPI API
#ifdef USEMPI
#include <mpi.h>
///Includes external global variables used for MPI version of code
#include "mpivar.h"
#endif

extern StrucLevelData *psldata;

#endif
