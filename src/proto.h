/*! \file proto.h
 *  \brief this file contains all function prototypes of the code

 \section TOC Function Type Table of contents
 */

#include "allvars.h"

#include "fofalgo.h"
#include "stf-fitting.h"

#ifndef STFPROTO_H
#define STFPROTO_H

//-- Prototypes

/// \name UI subroutines
/// see \ref ui.cxx for implementation
//@{

void usage(void);
void GetArgs(const int argc, char *argv[], Options &opt);
void GetParamFile(Options &opt);
void ConfigCheck(Options &opt);
void NOMASSCheck(Options &opt);

//@}

/// \name IO routines
/// see \ref io.cxx for implementation or for reading data see for example \ref gadgetio.cxx or \ref tipsyio.cxx
//@{

///simple check to see if file exists
bool FileExists(const char *fname);

//-- Read routines

///Reads the header information
Int_t ReadHeader(Options &opt);
///Reads particle data
void ReadData(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons=0);
///Read gadget file
void ReadGadget(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons=0);
///Read tipsy file
void ReadTipsy(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons=0);
#ifdef USEHDF
///Read HDF format
void ReadHDF(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons=0);
#endif
///Read ramses file
void ReadRamses(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons=0);
///Read Nchilada file
void ReadNchilada(Options &opt, vector<Particle> &Part, const Int_t nbodies,Particle *&Pbaryons, Int_t nbaryons=0);
///Adjust hydro particles/quantities to appropriate units
void AdjustHydroQuantities(Options &opt, vector<Particle> &Part, const Int_t nbodies);
///Adjust star particles/quantities to appropriate units
void AdjustStarQuantities(Options &opt, vector<Particle> &Part, const Int_t nbodies);
///Adjust BH particles/quantities to appropriate units
void AdjustBHQuantities(Options &opt, vector<Particle> &Part, const Int_t nbodies);

///Read local velocity density
void ReadLocalVelocityDensity(Options &opt, const Int_t nbodies, vector<Particle> &Part);
///Writes local velocity density of each particle to a file
void WriteLocalVelocityDensity(Options &opt, const Int_t nbodies, vector<Particle> &Part);


///Writes a tipsy formatted fof.grpfile
void WriteFOF(Options &opt, const Int_t nbodies, Int_t *pfof);
///Writes a pg list file (first in effective index order of input file(s), second is particle ids
void WritePGList(Options &opt, const Int_t ngroups, const Int_t ng, Int_t *numingroup, Int_t **pglist, Int_t *ids);
///Write catalog information (number of groups, number in groups, number of particles in groups, particle pids)
void WriteGroupCatalog(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, vector<Particle> &Part, Int_t nadditional=0);
///Write catalog information related to particle types relevant if different particle types are included in the grouping algorithm
void WriteGroupPartType(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, vector<Particle> &Part);
///Writes the bulk properties of the substructures
void WriteProperties(Options &opt, const Int_t ngroups, PropData *pdata);
///Writes the bulk properties of the substructures in a subfind like HDF5 format
void WriteSUBFINDProperties(Options &opt, const Int_t ngroups, PropData *pdata);
///Writes the structure hierarchy
//void WriteHierarchy(Options &opt, Int_t ngroups, int subflag=0);
void WriteHierarchy(Options &opt, const Int_t &ngroups, const Int_t &nhierarchy, const Int_t &nfield, Int_t *nsub, Int_t *parentgid, Int_t *stype,int subflag=0);

///Write Extended Output
void WriteExtendedOutput(Options &opt, Int_t numgroups, Int_t nbodies, PropData *pdata, vector<Particle> &p, Int_t * pfof);
#ifdef SWIFTINTERFACE
///Write catalog information related to the address of particles in groups in a swift snapshot.
void WriteSwiftExtendedOutput(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, vector<Particle> &Part);
#endif

///Write the configuartion options that the code used
void WriteVELOCIraptorConfig(Options &opt);
///Write the simulation info (which could use input files to overwrite passed configuration options)
void WriteSimulationInfo(Options &opt);
///Write the unit info
void WriteUnitInfo(Options &opt);

///Write particle ids of those within spherical overdensity of a field halo
void WriteSOCatalog(Options &opt, const Int_t ngroups, vector<Int_t> *SOpids, vector<int> *SOtypes=NULL);
///Write profiles
void WriteProfiles(Options &opt, const Int_t ngroups, PropData *pdata);
///Writes ROCKSTAR like output
//@{
//@}

///Writes AHF like output
//@{
///writes
void WritePGListIndex(Options &opt, const Int_t ngroups, const Int_t ng, Int_t *numingroup, Int_t **pglist);
///Writes the bulk properties of the substructures
void WriteAHFProperties(Options &opt, const Int_t ngroups, PropData *pdata);
//@}

///read a group file which can be useful if a halo search has already been run and one only wishes to
///identify substructures
//@{
Int_t ReadPFOF(Options &opt, Int_t nbodies, Int_t *pfof);
Int_t ReadFOFGroupBinary(Options &opt, Int_t nbodies, Int_t *pfof, Int_t *idtoindex, Int_t minid, Particle *p);
//@}

///Reads properties of a subvolume
void ReadCellValues(Options &opt, const Int_t nbodies, const Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp);
// ///Writes cell quantites
//void WriteCellValues(Options &opt, const Int_t nbodies, const Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp);
// ///Reads number of cells
//Int_t ReadCellNum(Options &opt);

///Print some info about cosmology
void PrintCosmology(Options &opt);
///Print some info about simulation
void PrintSimulationState(Options &opt);

#ifdef EXTENDEDHALOOUTPUT
//Write extended halo output for extraction of haloes from input files
void WriteExtendedOutput (Options &opt, Int_t numgroups, Int_t nbodies, PropData *pdata, Particle *p, Int_t * pfof);
#endif
//@}

/// \name Subroutines set up nonuniform grid of cells
/// see \ref bgfield.cxx for implementation
//@{

///Set up non-uniform grid structure using kd-tree
KDTree* InitializeTreeGrid(Options &opt, const Int_t nbodies, Particle *Part);
///Fill cells of grid from tree
void FillTreeGrid(Options &opt, const Int_t nbodies, const Int_t ngrid, KDTree *&tree, Particle *Part, GridCell* &grid);

//@}

/// \name Subroutines calculate mean velocity and velocity dispersion using grid
/// see \ref bgfield.cxx for implementation
//@{

///Calculate CM vel of cell
Coordinate *GetCellVel(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngrid, GridCell *grid);
///Calculate velocity dispersion tensor of cell
Matrix *GetCellVelDisp(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngrid, GridCell *grid, Coordinate *gvel);
//@}

/// \name Subroutines to calculate local velocity density
/// see \ref localfield.cxx for implementation
//@{

///Calculate local velocity density
void GetVelocityDensity(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree=NULL);
///sub interfaces depending on type of velocity density desired.
void GetVelocityDensityOld(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree);
///Velocity density where only particles in a halo (which is localised to an mpi domain) are within the tree
void GetVelocityDensityHaloOnlyDen(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree);
///exact velocity density, finds for each particle nearest physical neighbours and estimates velocity
void GetVelocityDensityExact(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree);
///optimised search for cosmological simulations
void GetVelocityDensityApproximative(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree);
//@}

/// \name Suboutines that compare local velocity density to background
/// see \ref localbgcomp.cxx for implementation
//@{

///Calculate logarithmic contrast between local velocity density and background velocity density
void GetDenVRatio(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp);
///Characterize the ratio distribution to estimate intrinsic scatter in this estimator
void DetermineDenVRatioDistribution(Options &opt,const Int_t nbodies, Particle *Part, Double_t &meanr,Double_t &sdlow,Double_t &sdhigh, int subleve=0);
///Calculate normalized residual
Int_t GetOutliersValues(Options &opt, const Int_t nbodies, Particle *Part, int sublevel=0);

//@}

/*! \name Search routines
    \brief see \ref search.cxx for implementation

    The search routines build a kd-tree and then search the tree to find particles that meet the fof criterion. \n

    If \b MPI is used then a tree is build for the local particle set and these particles are searched. Then
    \arg Once that is done, the node structure of each thread (with node upper most boundaries) is broadcast to all other processors
    so that a local thread can determine if a local particle can be linked to particles on the other threads.
    \arg This mpi_node structure is used to create a export list of the local particles that need to be broadcast to another thread \e j
    \arg After creating an export list of particles, these particles are sent and received by and to the appropriate threads.
    \arg Now the local thread searches it local particle set using the export particle positions to see if any links can be found.
    \n
    All these above steps must be repeated till no new links are found.
*/
//@{

///Search full system without finding outliers first
Int_t *SearchFullSet(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t &numgroups);
///Search the outliers
Int_t *SearchSubset(Options &opt, const Int_t nbodies, const Int_t nsubset, Particle *Partsubset, Int_t &numgroups, Int_t sublevel=0, Int_t *pnumcores=NULL);
///Search for subsubstructures
void SearchSubSub(Options &opt, const Int_t nsubset, vector<Particle> &Partsubset, Int_t *&pfof, Int_t &ngroup, Int_t &nhalos, PropData *pdata=NULL);
///Given a set of tagged core particles, assign surroundings
void HaloCoreGrowth(Options &opt, const Int_t nsubset, Particle *&Partsubset, Int_t *&pfof, Int_t *&pfofbg, Int_t &numgroupsbg, Double_t param[], vector<Double_t> &dispfac,
    int numactiveloops, vector<int> &corelevel, int nthreads);
///merge substructures and cores if phase-space positions overlap
void MergeSubstructuresCoresPhase(Options &opt, const Int_t nsubset, Particle *&Partsubset, Int_t *&pfof, Int_t &numsubs, Int_t &numcores);
///merge substructures if phase-space positions overlap
void MergeSubstructuresPhase(Options &opt, const Int_t nsubset, Particle *&Partsubset, Int_t *&pfof, Int_t &numgroups, Int_t &numsubs, Int_t &numcores);
///remove spurious dynamical substructures that comprise most of their host halo
void RemoveSpuriousDynamicalSubstructures(Options &opt, const Int_t nsubset, Int_t *&pfof, Int_t &numgroups, Int_t &numsubs, Int_t &numcores);
///Check significance of each group
int CheckSignificance(Options &opt, const Int_t nsubset, Particle *Partsubset, Int_t &numgroups, Int_t *numingroups, Int_t *pfof, Int_t **pglist);
///Search for Baryonic structures associated with dark matter structures in phase-space
Int_t* SearchBaryons(Options &opt, Int_t &nbaryons, Particle *&Pbaryons, const Int_t ndark, vector<Particle> &Partsubset, Int_t *&pfofdark, Int_t &ngroupdark, Int_t &nhalos, int ihaloflag=0, int iinclusive=0, PropData *phalos=NULL);
///Get the hierarchy of structures found
Int_t GetHierarchy(Options &opt, Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t *stype);
///Copy hierarchy to PropData structure
void CopyHierarchy(Options &opt, PropData *pdata,Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t *stype);
///Adjust particles in structures for period such that all substructure searches no longer have to run periodic searches
void AdjustStructureForPeriod(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t numgroups, Int_t *pfof);
///Adjust halo (object) positions back to inside the periodic volume.
void AdjustHaloPositionForPeriod(Options &opt, Int_t ngroup, Int_t *&numingroup, PropData *&pdata);
//@}

/// \name Extra routines used in iterative search
//@{

///for halo mergers check, get mean velocity dispersion
inline Double_t GetVelDisp(Particle *Part, Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///Create array for each group listing all previously unlinked particles that should now be linked.
inline void DetermineNewLinks(Int_t nsubset, Particle *Partsubset, Int_t *pfof, Int_t numgroups, Int_t &newlinks,
    Int_t *newlinksIndex, Int_t *numgrouplinksIndex, Int_t *nnID, Int_t ***newIndex);
///Create array for each group listing all possibly intergroup links
inline void DetermineGroupLinks(Int_t nsubset, Particle *Partsubset, Int_t *pfof, Int_t numgroups, Int_t &newlinks,
    Int_t *newlinksIndex, Int_t *numgrouplinksIndex, Int_t *nnID, Int_t ***newIndex);
///once Group links are found, reduce the data so that one can determine for each group, number of other groups linked, their gids and the number of that group linked
///these are stored in numgrouplinksIndex, intergroupgidIndex, & newintergroupIndex
inline void DetermineGroupMergerConnections(Particle *Partsubset, Int_t numgroups, Int_t *pfof, int *ilflag,
    Int_t *numgrouplinksIndex, Int_t *intergrouplinksIndex, Int_t *nnID, Int_t ***newIndex, Int_t ***newintergroupIndex, Int_t ***intergroupgidIndex);
///used to search particle list using tree to tag particles based on a comparison function (here distance information is used)
inline void SearchForNewLinks(Int_t nsubset, KDTree *tree, Particle *Partsubset, Int_t *pfof, FOFcompfunc &fofcmp, Double_t *param,
    Int_t newlinks, Int_t *newlinksIndex, Int_t **nnID, Double_t **dist2, int nthreads);
///used to search particle list using tree to tag particles based on a comparison function (here distance information is NOT used)
inline void SearchForNewLinks(Int_t nsubset, KDTree *tree, Particle *Partsubset, Int_t *pfof, FOFcompfunc &fofcmp, Double_t *param,
    Int_t newlinks, Int_t *newlinksIndex, Int_t **nnID, int nthreads);
///used to link new tag particles1
inline void LinkUntagged(Particle *Partsubset, Int_t numgroups, Int_t *pfof, Int_t *numingroup, Int_t **pglist,
    Int_t newlinks, Int_t *numgrouplinksIndex, Int_t **newIndex,
    Int_tree_t *Head, Int_tree_t *Next, Int_tree_t *GroupTail, Int_t *nnID);
///used to merge candidate substructure groups for iterative search
inline Int_t MergeGroups(Options &opt, Particle *Partsubset, Int_t numgroups, Int_t *pfof, Int_t *numingroup, Int_t *oldnumingroup,
    Int_t **pglist, Int_t *numgrouplinksIndex, Int_t ***intergroupgidIndex, Int_t ***newintergroupIndex, Int_t *intergrouplinksIndex,
    Int_tree_t *Head, Int_tree_t *Next, Int_tree_t *GroupTail, int *igflag, Int_t *nnID, Int_t &newlinks, Int_t *newlinksIndex);
///used to merge groups for halo merger check
inline Int_t MergeHaloGroups(Options &opt, Particle *Partsubset, Int_t numgroups, Int_t *pfof, Int_t *numingroup, Int_t *oldnumingroup,
    Int_t **pglist, Int_t *numgrouplinksIndex, Int_t ***intergroupgidIndex, Int_t ***newintergroupIndex, Int_t *intergrouplinksIndex,
    Int_tree_t *Head, Int_tree_t *Next, Int_tree_t *GroupTail, int *igflag, Int_t *nnID, Int_t &newlinks, Int_t *newlinksIndex);
//@}
/// \name For unbinding
/// see \ref unbind.cxx for implementation
//@{

///used for tree potential calculation (which is only used for large groups)
void GetNodeList(Node *np, Int_t &ncell, Node **nodelist, const Int_t bsize);
///used for tree walk in potential calculation
inline void MarkCell(Node *np, Int_t *marktreecell, Int_t *markleafcell, Int_t &ntreecell, Int_t &nleafcell, const Int_t bsize, Double_t *cR2max, Coordinate *cm, Double_t *cmtot, Coordinate xpos, Double_t eps2);

///Interface for unbinding proceedure
int CheckUnboundGroups(Options opt, const Int_t nbodies, Particle *Part, Int_t &ngroup, Int_t *&pfof, Int_t *numingroup=NULL, Int_t **pglist=NULL,int ireorder=1, Int_t *groupflag=NULL);
///check if group self-bound
int Unbind(Options &opt, Particle **gPartList, Int_t &numgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, int ireorder=1);
int Unbind(Options &opt, Particle *Part, Int_t &numgroups, Int_t *&numingroup, Int_t *&noffset, Int_t *&pfof);
///calculate the potential of an array of particles
void Potential(Options &opt, Int_t nbodies, Particle *Part, Double_t *potV);
void Potential(Options &opt, Int_t nbodies, Particle *Part);
void ParticleSubSample(Options &opt, const Int_t nbodies, Particle *&Part,
    Int_t &newnbodies, Particle *&newpart, double &mr);
void PotentialTree(Options &opt, Int_t nbodies, Particle *&Part, KDTree* &tree);
void PotentialInterpolate(Options &opt, const Int_t nbodies, Particle *&Part, Particle *&interolateparts, KDTree *&tree, double massratio, int nsearch);

void PotentialPP(Options &opt, Int_t nbodies, Particle *Part);
//@}

/// \name Routines to determine bulk quantities of halo and adjust halo
/// see \ref haloproperties.cxx for implementation
//@{

///Scales Linking lengths
void ScaleLinkingLengths(Options &opt,const Int_t nbodies, Particle *Part,Coordinate &cm,Coordinate &cmvel,Double_t Mtot);
///Adjust to cm of halo
void AdjusttoCM(const Int_t nbodies, Particle *Part, Coordinate &cm, Coordinate &cmvel, const Double_t Mtot, Double_t *rlim, Double_t &MaxVcirc, Double_t tol=1e-2);
///Calculate virial radius and related quantities
void GetVirialQuantities(const Int_t nbodies, Particle *Part, const Double_t mtot, Double_t *rlim, const Double_t rhoc, const Double_t virlevel, Double_t &Rvir, Double_t &Mvir, Int_t nenc, Double_t *Menc, Double_t *Renc);

//@}

/// \name Routines for cosmology related calculations
//@{
void CalcOmegak(Options &opt);
void CalcCriticalDensity(Options &opt, Double_t a);
void CalcBackgroundDensity(Options &opt, Double_t a);
void CalcVirBN98(Options &opt, Double_t a);
void CalcCosmoParams(Options &opt, Double_t a);
Double_t CalcGravitationalConstant(Options &opt);
Double_t CalcHubbleUnit(Options &opt);
Double_t GetHubble(Options &opt, Double_t a);
double GetInvaH(double a, void * params);
Double_t CalcCosmicTime(Options &opt, Double_t a1, Double_t a2);
//@}
/// \name Routines to calculate substructure properties and sort particles in a substructure according to some property
/// see \ref substructureproperties.cxx for implementation
//@{
///Get Centre of mass quantities
void GetCM(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
///Get the properties of the substructures and output the results
void GetProperties(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
///Get CM properties
//void GetCMProp(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
///Get max dist to several reference frames
void GetMaximumSizes(Options &opt, Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
///Adjust positions of bulk properties to desired reference
void AdjustHaloPositionRelativeToReferenceFrame(Options &opt, Int_t ngroup, Int_t *&numingroup, PropData *&pdata);
///Get NFW R200crit concentrations assuming an NFW profile and using Vmax/Rmax
void GetNFWConcentrations(Options &opt, Int_t ngroup, Int_t *&numingroup, PropData *&pdata);
///Get inclusive masses for field objects
void GetInclusiveMasses(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
///Get FOF masses for field objects
void GetFOFMass(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
/// Calculate FOF mass looping over groups once substructure search and have calculated properties
void GetFOFMass(Options &opt, Int_t ngroup, Int_t *&numingroup, PropData *&pdata);
///Get Spherical Overdensity Masses for field objects only. Assumes want to keep input order.
void GetSOMasses(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&numingroup, PropData *&pdata);
///simple routine to copy over mass information (useful for storing inclusive info)
void CopyMasses(Options &opt, const Int_t nhalos, PropData *&pold, PropData *&pnew);
///simple routine to reorder mass information based on number of particles when new remaining number of haloes < old halos
void ReorderInclusiveMasses(const Int_t &nold, const Int_t &nnew, Int_t *&numingroup, PropData *&pdata);
///Get Binding Energy
void GetBindingEnergy(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);

///Get Morphology properties (since this is for a particular system just use pointer interface)
void GetGlobalSpatialMorphology(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, Matrix& eigenvec, int imflag=0, int itype=-1, int iiterate=1);
///Calculate inertia tensor and eigvector
void CalcITensor(const Int_t n, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix &I, int itype);
///Calculate position dispersion tensor and eigvector
void CalcPosSigmaTensor(const Int_t n, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix &I, int itype=-1);
///Calculate velocity dispersion tensor and eigvector
void CalcVelSigmaTensor(const Int_t n, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix &I, int itype=-1);
///Calculate phase-space dispersion tensor and eigvector
void CalcPhaseSigmaTensor(const Int_t n, Particle *p, GMatrix &eigenvalues, GMatrix& eigenvec, GMatrix &I, int itype=-1);
///Calculate phase-space dispersion tensor
void CalcPhaseSigmaTensor(const Int_t n, Particle *p, GMatrix &I, int itype=-1);
///Calculate the reduced weighted inertia tensor used to determine the spatial morphology
void CalcMTensor(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p, int itype);
///Same as \ref CalcMTensor but include mass
void CalcMTensorWithMass(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p, int itype);
///Rotate particles to some coordinate frame
void RotParticles(const Int_t n, Particle *p, Matrix &R);
///get phase-space center-of-mass
GMatrix CalcPhaseCM(const Int_t n, Particle *p, int itype=-1);

///get concentration routines associted with finding concentrations via root finding
void CalcConcentration(PropData &p);
//root finding using half masses
double CalcConcentrationRootFindingRhalf(double , double);
//root finding using Vmax
double CalcConcentrationRootFindingVmax(double , double);
///wrappers for root finding used to get concentration
double mycNFW(double c, void *params);
///wrappers for root finding used to get concentration
double mycNFW_deriv(double c, void *params);
///wrappers for root finding used to get concentration
double mycNFW_fdf(double c, void *params, double*y,double *dy);
///wrappers for root finding used to get concentration
double mycNFWRhalf(double c, void *params);
///Calculate aperture quantities
void CalculateApertureQuantities(Options &opt, Int_t &ning, Particle *Part, PropData &pdata);
///determine the radial bin for calculating profiles
int GetRadialBin(Options &opt, Double_t rc, int &ibin);
///add a particle's properties to the appropriate radial bin.
void AddParticleToRadialBin(Options &opt, Particle *Pval, Double_t irnorm, int &ibin, PropData &pdata);
///add data to the appropriate radial bin
void AddDataToRadialBin(Options &opt, Double_t rval, Double_t massval,
#if defined(GASON) || defined(STARON) || defined(BHON)
    Double_t srfval, int typeval,
#endif
    Double_t irnorm, int &ibin, PropData &pdata);
void AddParticleToRadialBinInclusive(Options &opt, Particle *Pval, Double_t irnorm, int &ibin, PropData &pdata);
///add data to the appropriate radial bin
void AddDataToRadialBinInclusive(Options &opt, Double_t rval, Double_t massval,
#if defined(GASON) || defined(STARON) || defined(BHON)
    Double_t srfval, int typeval,
#endif
    Double_t irnorm, int &ibin, PropData &pdata);

///calculate extra hydro properties
void GetExtraHydroProperties(Options &opt, PropData &pdata, Int_t n, Particle *Pval);
///calculate extra star properties
void GetExtraStarProperties(Options &opt, PropData &pdata, Int_t n, Particle *Pval);
///calculate extra bh properties
void GetExtraBHProperties(Options &opt, PropData &pdata, Int_t n, Particle *Pval);
///calculate extra dm properties
void GetExtraDMProperties(Options &opt, PropData &pdata, Int_t n, Particle *Pval);

///calculate spherical overdensity from vector of radii, masses and indices
Int_t CalculateSphericalOverdensity(Options &opt, PropData &pdata,
    vector<Double_t> &radii, vector<Double_t> &masses, vector<Int_t> &indices,
    Double_t &m200val, Double_t &m200mval, Double_t &mBN98val, Double_t &virval, Double_t &m500val,
    vector<Double_t> &SOlgrhovals);
Int_t CalculateSphericalOverdensity(Options &opt, PropData &pdata,
    Int_t &numingroup, Particle *Part,
    Double_t &m200val, Double_t &m200mval, Double_t &mBN98val, Double_t &virval, Double_t &m500val,
    vector<Double_t> &SOlgrhovals);
void CalculateSphericalOverdensitySubhalo(Options &opt, PropData &pdata,
    Int_t &numingroup, Particle *Part,
    Double_t &m200val, Double_t &m200mval, Double_t &mBN98val, Double_t &virval, Double_t &m500val,
    vector<Double_t> &SOlgrhovals);
void CalculateSphericalOverdensityExclusive(Options &opt, PropData &pdata,
    Int_t &numingroup, Particle *Part,
    Double_t &m200val, Double_t &m200mval, Double_t &mBN98val, Double_t &virval, Double_t &m500val,
    vector<Double_t> &SOlgrhovals);

void SetSphericalOverdensityMasstoFlagValue(Options &opt, PropData &pdata);
void SetSphericalOverdensityMasstoTotalMass(Options &opt, PropData &pdata);
void SetSphericalOverdensityMasstoTotalMassExclusive(Options &opt, PropData &pdata);
///used to sort a pglist based on substructure binding energy
Int_t **SortAccordingtoBindingEnergy(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *numingroup, PropData *pdata, Int_t ioffset=0);
///used to calculate properties and ignores keeping particle order, assumes particle PID information meaningless
void CalculateHaloProperties(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *numingroup, PropData *pdata);
///Get number of substructures in a given structure
Int_t *GetSubstructureNum(Int_t ngroups);
///Get parent structure id of substructures
Int_t *GetParentID(Int_t ngroups);
//@}

#ifdef USEOPENMP
/// \name OpenMP Search routines
/// see \ref omproutines.cxx for implementation
//@{
///build the openmp domains
OMP_Domain *OpenMPBuildDomains(Options &opt, const Int_t numompregions, KDTree *&tree, const Double_t rdist);

///build the local OpenMP trees that are used to search the domain
KDTree **OpenMPBuildLocalTrees(Options &opt, const Int_t numompregions, vector<Particle> &Part, OMP_Domain *ompdomain, Double_t *period);

///returns whether overlap is found
int OpenMPSearchForOverlap(Double_t xsearch[3][2], Double_t bnd[3][2], Double_t period);

///returns whether overlap is found
int OpenMPSearchForOverlap(Particle &Part, Double_t bnd[3][2], Double_t rdist, Double_t period);

///returns whether region fully contained within domain
int OpenMPInDomain(Double_t xsearch[3][2], Double_t bnd[3][2]);

///returns whether particles search region fully contained within domain
int OpenMPInDomain(Particle &Part, Double_t bnd[3][2], Double_t rdist);

///Search all OpenMP domains using FOF algorithm
Int_t OpenMPLocalSearch(Options &opt,
    const Int_t nbodies, vector<Particle> &Part, Int_t * &pfof, Int_t *&storeorgIndex,
    Int_tree_t *&Head, Int_tree_t *&Next,
    KDTree **&tree3dfofomp, Double_t *param, const Double_t rdist, const Int_t ompminsize, FOFcompfunc fofcomp,
    const Int_t numompregions, OMP_Domain *&ompdomain);

///determine particle to import from other OpenMP domains
OMP_ImportInfo *OpenMPImportParticles(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t * &pfof, Int_t *&storetype,
    const Int_t numompregions, OMP_Domain *&ompdomain, const Double_t rdist,
    Int_t *&omp_nrecv_total, Int_t *&omp_nrecv_offset, Int_t &omp_import_total);

///link across mpi domains
void OpenMPLinkAcross(Options &opt,
    Int_t nbodies, vector<Particle> &Part, Int_t * &pfof,
    Int_t *&storetype, Int_tree_t *&Head, Int_tree_t *&Next,
    Double_t *param, FOFcheckfunc &fofcheck,
    const Int_t numompregions, OMP_Domain *&ompdomain, KDTree **tree3dfofomp,
    Int_t *&omp_nrecv_total, Int_t *&omp_nrecv_offset, OMP_ImportInfo* &ompimport);

///resorts particles and group id values after OpenMP search
Int_t OpenMPResortParticleandGroups(Int_t nbodies, vector<Particle> &Part, Int_t *&pfof, Int_t minsize);

///sets the head/next arrays based on the current particle order and the current pfof array
void OpenMPHeadNextUpdate(const Int_t nbodies, vector<Particle> &Part, const Int_t numgroups, Int_t *&pfof, Int_tree_t *&Head, Int_tree_t *&Next);
#endif

#ifdef USEMPI
/// \name MPI Domain Decomposition routines
/// see \ref mpiroutines.cxx for implementation
/// for format specific stuff see related routines like \ref mpigadgetio.cxx or \ref mpitipsyio.cxx
//@{

///Update config options across mpi tasks of particles to load
void MPIUpdateUseParticleTypes(Options &opt);

///Domain decomposition of system
void MPIInitialDomainDecomposition(Options &opt);
///Domain extent
void MPIDomainExtent(Options &opt);
///domain decomposition
void MPIDomainDecomposition(Options &opt);
///z-curve based mesh decomposition
void MPIInitialDomainDecompositionWithMesh(Options &opt);
///z-curve repartitioning of cells
bool MPIRepartitionDomainDecompositionWithMesh(Options &opt);

///Determine Domain Extent for tipsy input
void MPIDomainExtentTipsy(Options &opt);
///Determine Domain for tipsy input
void MPIDomainDecompositionTipsy(Options &opt);

///Determine Domain Extent for Gadget input
void MPIDomainExtentGadget(Options &opt);
///Determine Domain for Gadget input
void MPIDomainDecompositionGadget(Options &opt);

///Determine Domain Extent for Ramses input
void MPIDomainExtentRAMSES(Options &opt);
///Determine Domain for Ramses input
void MPIDomainDecompositionRAMSES(Options &opt);

#ifdef USEHDF
///Determine Domain Extent for HDF input
void MPIDomainExtentHDF(Options &opt);
///Determine Domain for Gadget input
void MPIDomainDecompositionHDF(Options &opt);
#endif

///Determine Domain Extent for Nchilada input
void MPIDomainExtentNchilada(Options &opt);
///Determine Domain for Gadget input
void MPIDomainDecompositionNchilada(Options &opt);

///determine what processor a particle is sent to based on domain decomposition
int MPIGetParticlesProcessor(Options &opt, const Double_t,const Double_t,const Double_t);
/// Determine number of local particles wrapper
void MPINumInDomain(Options &opt);

///determine which mpi processes read input files
void MPIDistributeReadTasks(Options&opt, int *&ireadtask, int*&readtaskID);
///set which file a given task will read
int MPISetFilesRead(Options&opt, int *&ireadfile, int *&ireadtask);

///generic init of write communicator to mpi world;
void MPIInitWriteComm();
///determine how to group mpi threads together when writing
void MPIBuildWriteComm(Options &opt);
///free communicator if needed
void MPIFreeWriteComm();

/// Determine number of local particles for tipsy
void MPINumInDomainTipsy(Options &opt);
/// Determine number of local particles for gadget
void MPINumInDomainGadget(Options &opt);
/// Determine number of local particles for ramses
void MPINumInDomainRAMSES(Options &opt);
#ifdef USEHDF
/// Determine number of local particles for HDF
void MPINumInDomainHDF(Options &opt);
#endif
/// Determine number of local particles for Nchilada
void MPINumInDomainNchilada(Options &opt);

///adjust the domain boundaries to code units
void MPIAdjustDomain(Options &opt);
///adjust the domain boundaries to code units
void MPIAdjustDomainWithMesh(Options &opt);
///determine if the search domain of a particle overlaps another mpi domain
int MPISearchForOverlap(Particle &Part, Double_t &rdist);
///determine if the search domain of a particle overlaps another mpi domain
int MPISearchForOverlap(Coordinate &x, Double_t &rdist);
///determine if the search domain of a particle overlaps another mpi domain using the SWIFT mesh
int MPISearchForOverlapUsingMesh(Options &opt, Particle &Part, Double_t &rdist);
///determine if the search domain of a coordinate overlaps another mpi domain using the SWIFT mesh
int MPISearchForOverlapUsingMesh(Options &opt, Coordinate &x, Double_t &rdist);
///determine if the search domain overlaps another mpi domain
int MPISearchForOverlap(Double_t xsearch[3][2]);
int MPISearchForOverlapUsingMesh(Options &opt, Double_t xsearch[3][2]);

///determine if search domain overlaps domain
int MPIInDomain(Double_t xsearch[3][2], Double_t bnd[3][2]);

///determine list of cells of a mesh within a search domain
vector<int> MPIGetCellListInSearchUsingMesh(Options &opt, Double_t xsearch[3][2], bool ignorelocalcells=true);
vector<int> MPIGetCellNodeIDListInSearchUsingMesh(Options &opt, Double_t xsearch[3][2]);
//@}

/// \name MPI send/recv related routines when reading input data
/// see \ref mpiroutines.cxx for implementation
//@{

///adds particles to appropriate send buffers and initiates sends if necessary.
void MPIAddParticletoAppropriateBuffer(Options &opt, const int &ibuf, Int_t ibufindex, int *&ireadtask, const Int_t &Bufsize, Int_t *&Nbuf, Particle *&Pbuf, Int_t &numpart, Particle *Part, Int_t *&Nreadbuf, vector<Particle>*&Preadbuf);
///Send particle information from read threads to non read threads using MPI_COMM_WORLD
void MPISendParticlesFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID);
///recv particle data from read threads
void MPIReceiveParticlesFromReadThreads(Options &opt, Particle *&Pbuf, Particle *Part, int *&readtaskID, int *&irecv, int *&mpi_irecvflag, Int_t *&Nlocalthreadbuf, MPI_Request *&mpi_request, Particle *&Pbaryons);
///Send/recv particle data read from input files between the various read threads;
void MPISendParticlesBetweenReadThreads(Options &opt, Particle *&Pbuf, Particle *Part, Int_t *&nreadoffset, int *&ireadtask, int *&readtaskID, Particle *&Pbaryons, Int_t *&mpi_nsend_baryon);
///Send/recv particle data stored in vector using the read thread communication domain
void MPISendParticlesBetweenReadThreads(Options &opt, vector<Particle> *&Pbuf, Particle *Part, int *&ireadtask, int *&readtaskID, Particle *&Pbaryons, MPI_Comm &mpi_read_comm, Int_t *&mpi_nsend_readthread, Int_t *&mpi_nsend_readthread_baryon);

///Interrupt send of particle information to destination taskID using MPI_COMM_WORLD
void MPIISendParticleInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID, int tag, MPI_Request &rqst);
///Receive Particle information send with specific tag
void MPIReceiveParticleInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sendingTaskID, int tag);

///Send hydro information from read threads to non read threads using MPI_COMM_WORLD
void MPISendHydroInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID);
///Send star information from read threads to non read threads using MPI_COMM_WORLD
void MPISendStarInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID);
///Send bh information from read threads to non read threads using MPI_COMM_WORLD
void MPISendBHInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID);
///Send extra dm information from read threads to non read threads using MPI_COMM_WORLD
void MPISendExtraDMInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID);

///Receive hydro information from read threads using MPI_COMM_WORLD
void MPIReceiveHydroInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID);
///Receive star information from read threads using MPI_COMM_WORLD
void MPIReceiveStarInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID);
///Receive bh information from read threads using MPI_COMM_WORLD
void MPIReceiveBHInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID);
///Receive extra dm information from read threads using MPI_COMM_WORLD
void MPIReceiveExtraDMInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID);

///Interrupt send of hydro information to destination taskID using MPI_COMM_WORLD
void MPIISendHydroInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID, int tag, MPI_Request &rqst);
///Interrupt send of star information to destination taskID using MPI_COMM_WORLD
void MPIISendStarInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID, int tag, MPI_Request &rqst);
///Interrupt send of BH information to destination taskID using MPI_COMM_WORLD
void MPIISendBHInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID, int tag, MPI_Request &rqst);
///Interrupt send of extra dm information to destination taskID using MPI_COMM_WORLD
void MPIISendExtraDMInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID, int tag, MPI_Request &rqst);

///Receive Hydro information send with specific tag
void MPIReceiveHydroInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sendingTaskID, int tag);
///Receive star information send with specific tag
void MPIReceiveStarInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sendingTaskID, int tag);
///Receive BH information send with specific tag
void MPIReceiveBHInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sendingTaskID, int tag);
///Receive Extra DM information send with specific tag
void MPIReceiveExtraDMInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sendingTaskID, int tag);

///Send/Receive hydro information between read threads using the MPI communicator
void MPISendReceiveHydroInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive star information between read threads using the MPI communicator
void MPISendReceiveStarInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive BH information between read threads using the MPI communicator
void MPISendReceiveBHInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive BH information between read threads using the MPI communicator
void MPISendReceiveExtraDMInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm);

///strip off extra information stored in unique pointers if not necessary to send info
void MPIStripExportParticleOfExtraInfo(Options &opt, Int_t n, Particle *Part);

///Filling extra buffers with hydro data for particles that are to be exported
void MPIFillBuffWithHydroInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff=false);
///Filling extra buffers with star data for particles that are to be exported
void MPIFillBuffWithStarInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff=false);
///Filling extra buffers with bh data for particles that are to be exported
void MPIFillBuffWithBHInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff=false);
///Filling extra buffers with extra dm data for particles that are to be exported
void MPIFillBuffWithExtraDMInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff=false);

///Send/Receive hydro information between read threads using the MPI communicator
void MPISendReceiveBuffWithHydroInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive star information between read threads using the MPI communicator
void MPISendReceiveBuffWithStarInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive bh information between read threads using the MPI communicator
void MPISendReceiveBuffWithBHInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive extra dm information between read threads using the MPI communicator
void MPISendReceiveBuffWithExtraDMInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);

///Filling extra buffers with hydro data for particles that are to be exported as part of
///a FOF group.
void MPIFillFOFBuffWithHydroInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport = true);
///Filling extra buffers with hydro data for particles that are to be exported as part of
///a Star group.
void MPIFillFOFBuffWithStarInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport = true);
///Filling extra buffers with hydro data for particles that are to be exported as part of
///a BH group.
void MPIFillFOFBuffWithBHInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport = true);
///Filling extra buffers with hydro data for particles that are to be exported as part of
///a Extra DM group.
void MPIFillFOFBuffWithExtraDMInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport = true);

///Send/Receive hydro information between read threads using the MPI communicator
///Using a FOF filled buffer
void MPISendReceiveFOFHydroInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive star information between read threads using the MPI communicator
void MPISendReceiveFOFStarInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive bh information between read threads using the MPI communicator
void MPISendReceiveFOFBHInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
///Send/Receive extra dm information between read threads using the MPI communicator
void MPISendReceiveFOFExtraDMInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indices, vector<float> &propbuff, int recvTask, int tag, MPI_Comm &mpi_comm);
//@}

/// \name MPI search related routines
/// see \ref mpiroutines.cxx for implementation
//@{

///Set the mpi task ID of a particle to the local mpi thread or task number.
short_mpi_t *MPISetTaskID(const Int_t nbodies);
/// Adjust local group ids so that mpi threads are all offset from one another unless a particle belongs to group zero (ie: completely unlinked)
void MPIAdjustLocalGroupIDs(const Int_t nbodies, Int_t *pfof);
///Determine number of particles that need to be exported to another mpi thread from local mpi thread based on rdist
void MPIGetExportNum(const Int_t nbodies, Particle *Part, Double_t rdist);
///Determine number of particles that need to be exported to another mpi thread from local mpi thread based on rdist using the SWIFT mesh
void MPIGetExportNumUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t rdist);
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on rdist
void MPIBuildParticleExportList(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Double_t rdist);
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on rdist using the SWIFT mesh
void MPIBuildParticleExportListUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Double_t rdist);
///Link groups across MPI threads using a physical search
Int_t MPILinkAcross(const Int_t nbodies, KDTree *&tree, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Int_tree_t *&Head, Int_tree_t *&Next, Double_t rdist2);
///Link groups across MPI threads using criterion
Int_t MPILinkAcross(const Int_t nbodies, KDTree *&tree, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Int_tree_t *&Head, Int_tree_t *&Next, Double_t rdist2, FOFcompfunc &cmp, Double_t *params);
///Link groups across MPI threads checking particle types
Int_t MPILinkAcross(const Int_t nbodies, KDTree *&tree, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Int_tree_t *&Head, Int_tree_t *&Next, Double_t rdist2, FOFcheckfunc &check, Double_t *params);
///update export list after after linking across
void MPIUpdateExportList(const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len);
///localize groups to a single mpi thread
Int_t MPIGroupExchange(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof);
///Determine the local number of groups and their sizes (groups must be local to an mpi thread)
Int_t MPICompileGroups(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t minsize);
///similar to \ref MPIGroupExchange but optimised for separate baryon search, assumes only looking at baryons
Int_t MPIBaryonGroupExchange(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof);
///similar to \ref MPICompileGroups but optimised for separate baryon search, assumes only looking at baryons
Int_t MPIBaryonCompileGroups(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t minsize, int iorder=1);
///localize baryons particle members of groups to a single mpi thread
///Collect FOF from all
void MPICollectFOF(const Int_t nbodies, Int_t *&pfof);
///comparison function to order particles for export
int fof_export_cmp(const void *a, const void *b);
///comparison function to order particles for export and fof group localization.
int fof_id_cmp(const void *a, const void *b);
///comparison function to order particles for export
bool fof_export_cmp_vec(const fofdata_in &a, const fofdata_in &b);
///comparison function to order particles for export and fof group localization.
bool fof_id_cmp_vec(const fofid_in &a, const fofid_in &b);
///similar to \ref MPIBuildParticleExportList but specific interface for baryon search
void MPIBuildParticleExportBaryonSearchList(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t *ids, Int_t *numingroup, Double_t rdist);
///search local baryons with exported particle list.
Int_t MPISearchBaryons(const Int_t nbaryons, Particle *&Pbaryons, Int_t *&pfofbaryons, Int_t *numingroup, Double_t *localdist, Int_t nsearch, Double_t *param, Double_t *period);
///localize the baryons to the mpi thread on which their associated DM group exists.
Int_t MPIBaryonExchange(const Int_t nbaryons, Particle *&Pbaryons, Int_t *&pfofbaryons);
//@}

/// \name MPI NN determination related routines
/// see \ref mpiroutines.cxx for implementation
//@{

///Effectively is an allgather for the grid data so that particles can find nearest cells and use appropriate nearest neighbouring cells for calculating estimated background velocity density function
void MPIBuildGridData(const Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp);
///Determine number of particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIGetNNExportNum(const Int_t nbodies, Particle *Part, Double_t *rdist);
///Determine number of particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIGetNNExportNumUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t *rdist);
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIBuildParticleNNExportList(const Int_t nbodies, Particle *Part, Double_t *rdist);
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIBuildParticleNNExportListUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t *rdist);
///Determine number of local particles that need to be exported back based on ball search.
void MPIGetNNImportNum(const Int_t nbodies, KDTree *tree, Particle *Part, int iallflag=1);
///Determine local particles that need to be exported back based on ball search.
Int_t MPIBuildParticleNNImportList(Options &opt, const Int_t nbodies, KDTree *tree, Particle *Part, int iallflag=1, bool iSOcalc = false);
///comparison function to order particles for export
int nn_export_cmp(const void *a, const void *b);
///comparison function to order particles for export
bool nn_export_cmp_vec(const nndata_in &a, const nndata_in &b);
///Determine number of halos whose search regions overlap other mpi domains
vector<bool> MPIGetHaloSearchExportNum(const Int_t ngroups, PropData *&pdata, vector<Double_t> &rdist);
///Determine number of halos whose search regions overlap other mpi domains using swift mesh mpi decomposition
vector<bool> MPIGetHaloSearchExportNumUsingMesh(Options &opt, const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist);
///Build the export list of halo positions and search distances
void MPIBuildHaloSearchExportList(const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist, vector<bool> &halooverlap);
///Build the export list of halo positions and search distances using swift mesh mpi decomposition
void MPIBuildHaloSearchExportListUsingMesh(Options &opt, const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist, vector<bool> &halooverlap);
///Determine number of imported particles based on halo search regions
void MPIGetHaloSearchImportNum(const Int_t nbodies, KDTree *tree, Particle *Part);
///Builds the import list of particles based on halo positions
Int_t MPIBuildHaloSearchImportList(Options &opt, const Int_t nbodies, KDTree *tree, Particle *Part);
#ifdef SWIFTINTERFACE
///Exchange Particles so that particles in group are back original swift task
void MPISwiftExchange(vector<Particle> &Part);
#endif
//@}
#endif

/// \name Extra routines for building arrays that access the particle data or reorder arrays
/// see \ref buildandsortarrays.cxx for implementation
//@{

///build group size array
Int_t *BuildNumInGroup(const Int_t nbodies, const Int_t numgroups, Int_t *pfof);
///build group size array of particles of a specific type
Int_t *BuildNumInGroupTyped(const Int_t nbodies, const Int_t numgroups, Int_t *pfof, Particle *Part, int type);
///build array such that array is pglist[group][]={particle list}
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof);
///build array such that array is pglist[group][]={particle list} but only for particles of a specific type
///should be matched with \ref BuildNumInGroupTyped
Int_t **BuildPGListTyped(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *Part, int type);
///build the group particle index list (doesn't assume particles are in ID order and stores index of particle)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *Part);
///build pglist but doesn't assume particles are in ID order
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Int_t *ids);
///build the group particle arrays need for unbinding procedure
Particle **BuildPartList(const Int_t numgroups, Int_t *numingroup, Int_t **pglist, Particle* Part, bool ikeepextrainfo = false);
///build a particle list subset using array of indices
Particle *BuildPart(Int_t numingroup, Int_t *pglist, Particle* Part, bool ikeepextrainfo = false);
///build the Head array which points to the head of the group a particle belongs to
Int_tree_t *BuildHeadArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///build the Next array which points to the next particle in the group
Int_tree_t *BuildNextArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///build the Len array which stores the length of the group a particle belongs to
Int_tree_t *BuildLenArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///build the GroupTail array which stores the Tail of a group
Int_tree_t *BuildGroupTailArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///sort particles according to the group value (or technically any integer array) unique to each group and return an array of offsets to access the particle array via their group
Int_t *BuildNoffset(const Int_t nbodies, Particle *Part, Int_t numgroups,Int_t *numingroup, Int_t *sortval, Int_t ioffset=0);
///reorder groups from largest to smallest
void ReorderGroupIDs(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist);
///reorder groups from largest to smallest not assuming particles are in id order
void ReorderGroupIDs(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Particle *P);
///reorder groups by value
void ReorderGroupIDsbyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value);
///reorder groups and associated integer group data by value
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Int_t *gdata);
///reorder groups and associated double group data by value
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Double_t *gdata);
///reorder groups and associated integer group data by value
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Double_t *value, Int_t *gdata);
///reorder groups and associated double group data by value
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Double_t *value, Double_t *gdata);
///reorder groups and the associated property data by value
void ReorderGroupIDsAndHaloDatabyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, PropData *pdata);
//@}


/// \name Extra utility routines
/// see \ref utilities.cxx for implementation
//@{

int CompareInt(const void *, const void *);
///Get memory use
void GetMemUsage(Options &opt, string callingfunction, bool printreport);
///Get memory use
void GetMemUsage(string callingfunction, bool printreport);
///Init memory log
void InitMemUsageLog(Options &opt);
///get a time
std::chrono::time_point<std::chrono::high_resolution_clock> MyGetTime();
double MyElapsedTime(std::chrono::time_point<std::chrono::high_resolution_clock> before);
//@}

/// \name Compilation functions
/// Functions defined when certain compilation options enabled. Allows easy check
///of a library file produced. see \ref utilities.cxx for implementation
//@{
#ifdef NOMASS
extern "C" void VR_NOMASS();
#endif
#ifdef GASON
extern "C" void VR_GASON();
#endif
#ifdef STARON
extern "C" void VR_STARON();
#endif
#ifdef BHON
extern "C" void VR_BHON();
#endif
#ifdef USEMPI
extern "C" void VR_MPION();
#endif
#ifdef USEOPENMP
extern "C" void VR_OPENMPON();
#endif
#ifdef HIGHRES
extern "C" void VR_ZOOMSIMON();
#endif
#ifdef USEHDF
extern "C" void VR_HDFON();
#ifdef USEPARALLELHDF
extern "C" void VR_PARALLELHDFON();
#endif
#endif
//@}

#endif
