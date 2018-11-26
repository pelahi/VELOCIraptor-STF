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
inline void ConfigCheck(Options &opt);

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
///Writes the structure hierarchy
//void WriteHierarchy(Options &opt, Int_t ngroups, int subflag=0);
void WriteHierarchy(Options &opt, const Int_t &ngroups, const Int_t &nhierarchy, const Int_t &nfield, Int_t *nsub, Int_t *parentgid, Int_t *stype,int subflag=0);

///Write Extended Output
void WriteExtendedOutput(Options &opt, Int_t numgroups, Int_t nbodies, PropData *pdata, vector<Particle> &p, Int_t * pfof);

///Write the configuartion options that the code used
void WriteVELOCIraptorConfig(Options &opt);
///Write the simulation info (which could use input files to overwrite passed configuration options)
void WriteSimulationInfo(Options &opt);
///Write the unit info
void WriteUnitInfo(Options &opt);
///Write particle ids of those within spherical overdensity of a field halo
void WriteSOCatalog(Options &opt, const Int_t ngroups, vector<Int_t> *SOpids);
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
///Check significance of each group
int CheckSignificance(Options &opt, const Int_t nsubset, Particle *Partsubset, Int_t &numgroups, Int_t *numingroups, Int_t *pfof, Int_t **pglist);
///Search for Baryonic structures associated with dark matter structures in phase-space
Int_t* SearchBaryons(Options &opt, Int_t &nbaryons, Particle *&Pbaryons, const Int_t ndark, vector<Particle> &Partsubset, Int_t *&pfofdark, Int_t &ngroupdark, Int_t &nhalos, int ihaloflag=0, int iinclusive=0, PropData *phalos=NULL);
///Get the hierarchy of structures found
Int_t GetHierarchy(Options &opt, Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t *stype);
///Copy hierarchy to PropData structure
void CopyHierarchy(Options &opt, PropData *pdata,Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t *stype);
///Adjust Structure search for period such that all substructure searches no longer have to run periodic searches
void AdjustStructureForPeriod(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t numgroups, Int_t *pfof);
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

/// \name Routines to calculate substructure properties and sort particles in a substructure according to some property
/// see \ref substructureproperties.cxx for implementation
//@{
///Get the properties of the substructures and output the results
void GetProperties(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *numingroup=NULL, Int_t **pglist=NULL);
///Get CM properties
void GetCMProp(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
///Get inclusive masses for field objects
void GetInclusiveMasses(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *&numingroup, PropData *&pdata, Int_t *&noffset);
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
void GetConcentration(PropData &p);
///wrappers for root finding used to get concentration
double mycNFW(double c, void *params);
///wrappers for root finding used to get concentration
double mycNFW_deriv(double c, void *params);
///wrappers for root finding used to get concentration
double mycNFW_fdf(double c, void *params, double*y,double *dy);

///used to sort a pglist based on substructure binding energy
Int_t **SortAccordingtoBindingEnergy(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *numingroup, PropData *pdata, Int_t ioffset=0);
///used to calculate properties and ignores keeping particle order, assumes particle PID information meaningless
void CalculateHaloProperties(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *&pfof, Int_t *numingroup, PropData *pdata);
///Get number of substructures in a given structure
Int_t *GetSubstructureNum(Int_t ngroups);
///Get parent structure id of substructures
Int_t *GetParentID(Int_t ngroups);
//@}

#ifdef USEMPI
/// \name MPI Domain Decomposition routines
/// see \ref mpiroutines.cxx for implementation
/// for format specific stuff see related routines like \ref mpigadgetio.cxx or \ref mpitipsyio.cxx
//@{

///Domain decomposition of system
void MPIInitialDomainDecomposition();
///Domain extent
void MPIDomainExtent(Options &opt);
///domain decomposition
void MPIDomainDecomposition(Options &opt);

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
int MPIGetParticlesProcessor(const Double_t,const Double_t,const Double_t);
/// Determine number of local particles wrapper
void MPINumInDomain(Options &opt);

///determine which mpi processes read input files
void MPIDistributeReadTasks(Options&opt, int *&ireadtask, int*&readtaskID);
///set which file a given task will read
int MPISetFilesRead(Options&opt, int *&ireadfile, int *&ireadtask);

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
void MPIAdjustDomain(Options opt);
///determine if the search domain of a particle overlaps another mpi domain
int MPISearchForOverlap(Particle &Part, Double_t &rdist);
#ifdef SWIFTINTERFACE
///determine if the search domain of a particle overlaps another mpi domain using the SWIFT mesh
int MPISearchForOverlapUsingMesh(Options &opt, Particle &Part, Double_t &rdist);
#endif
///determine if the search domain overlaps another mpi domain
int MPISearchForOverlap(Double_t xsearch[3][2]);
///determine if search domain overlaps domain
int MPIInDomain(Double_t xsearch[3][2], Double_t bnd[3][2]);

///determine list of cells of a mesh within a search domain
vector<int> MPIGetCellListInSearchUsingMesh(Options &opt, Double_t xsearch[3][2], bool ignorelocalcells=true);
//@}

/// \name MPI send/recv related routines when reading input data
/// see \ref mpiroutines.cxx for implementation
//@{
///adds particles to appropriate send buffers and initiates sends if necessary.
void MPIAddParticletoAppropriateBuffer(const int &ibuf, Int_t ibufindex, int *&ireadtask, const Int_t &Bufsize, Int_t *&Nbuf, Particle *&Pbuf, Int_t &numpart, Particle *Part, Int_t *&Nreadbuf, vector<Particle>*&Preadbuf);
///recv particle data from read threads
void MPIReceiveParticlesFromReadThreads(Options &opt, Particle *&Pbuf, Particle *Part, int *&readtaskID, int *&irecv, int *&mpi_irecvflag, Int_t *&Nlocalthreadbuf, MPI_Request *&mpi_request, Particle *&Pbaryons);
///Send/recv particle data read from input files between the various read threads;
void MPISendParticlesBetweenReadThreads(Options &opt, Particle *&Pbuf, Particle *Part, Int_t *&nreadoffset, int *&ireadtask, int *&readtaskID, Particle *&Pbaryons, Int_t *&mpi_nsend_baryon);
///Send/recv particle data stored in vector using the read thread communication domain
void MPISendParticlesBetweenReadThreads(Options &opt, vector<Particle> *&Pbuf, Particle *Part, int *&ireadtask, int *&readtaskID, Particle *&Pbaryons, MPI_Comm &mpi_read_comm, Int_t *&mpi_nsend_readthread, Int_t *&mpi_nsend_readthread_baryon);
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
#ifdef SWIFTINTERFACE
///Determine number of particles that need to be exported to another mpi thread from local mpi thread based on rdist using the SWIFT mesh
void MPIGetExportNumUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t rdist);
#endif
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on rdist
void MPIBuildParticleExportList(const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Double_t rdist);
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
Int_t MPIGroupExchange(const Int_t nbodies, Particle *Part, Int_t *&pfof);
///Determine the local number of groups and their sizes (groups must be local to an mpi thread)
Int_t MPICompileGroups(const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t minsize);
///similar to \ref MPIGroupExchange but optimised for separate baryon search, assumes only looking at baryons
Int_t MPIBaryonGroupExchange(const Int_t nbodies, Particle *Part, Int_t *&pfof);
///similar to \ref MPICompileGroups but optimised for separate baryon search, assumes only looking at baryons
Int_t MPIBaryonCompileGroups(const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t minsize, int iorder=1);
///localize baryons particle members of groups to a single mpi thread
///Collect FOF from all
void MPICollectFOF(const Int_t nbodies, Int_t *&pfof);
///comparison function to order particles for export
int fof_export_cmp(const void *a, const void *b);
///comparison function to order particles for export and fof group localization.
int fof_id_cmp(const void *a, const void *b);
///similar to \ref MPIBuildParticleExportList but specific interface for baryon search
void MPIBuildParticleExportBaryonSearchList(const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t *ids, Int_t *numingroup, Double_t rdist);
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
#ifdef SWIFTINTERFACE
///Determine number of particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIGetNNExportNumUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t *rdist);
#endif
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIBuildParticleNNExportList(const Int_t nbodies, Particle *Part, Double_t *rdist);
#ifdef SWIFTINTERFACE
///Determine and send particles that need to be exported to another mpi thread from local mpi thread based on array of distances for each particle for NN search
void MPIBuildParticleNNExportListUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t *rdist);
#endif
///Determine number of local particles that need to be exported back based on ball search.
void MPIGetNNImportNum(const Int_t nbodies, KDTree *tree, Particle *Part);
///Determine local particles that need to be exported back based on ball search.
Int_t MPIBuildParticleNNImportList(const Int_t nbodies, KDTree *tree, Particle *Part, int iallflag=true);
///comparison function to order particles for export
int nn_export_cmp(const void *a, const void *b);
///Determine number of halos whose search regions overlap other mpi domains
vector<bool> MPIGetHaloSearchExportNum(const Int_t ngroups, PropData *&pdata, vector<Double_t> &rdist);
///Build the export list of halo positions and search distances
void MPIBuildHaloSearchExportList(const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist, vector<bool> &halooverlap);
///Determine number of imported particles based on halo search regions
void MPIGetHaloSearchImportNum(const Int_t nbodies, KDTree *tree, Particle *Part);
///Builds the import list of particles based on halo positions
Int_t MPIBuildHaloSearchImportList(const Int_t nbodies, KDTree *tree, Particle *Part);
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
Particle **BuildPartList(const Int_t numgroups, Int_t *numingroup, Int_t **pglist, Particle* Part);
///build a particle list subset using array of indices
Particle *BuildPart(Int_t numingroup, Int_t *pglist, Particle* Part);
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

///Get current time in Milliseconds
int GetMilliCount();
///Get span in milliseconds
int GetMillSpan(int );
int CompareInt(const void *, const void *);
///get a time
double MyGetTime();
//@}

#endif
