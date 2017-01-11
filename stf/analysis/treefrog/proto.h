/*! \file proto.h
 *  \brief this file contains all function prototypes of the code
 
 \section TOC Function Type Table of contents
 */

#include "allvars.h"

#ifndef HMTPROTO_H
#define HMTPROTO_H

//-- Prototypes

/// \name UI subroutines
/// see \ref ui.cxx for implementation
//@{

void usage(void);
void GetArgs(const int argc, char *argv[], Options &opt);
inline void ConfigCheck(Options &opt);

//@}

/// \name IO routines
/// see \ref io.cxx for implementation
//@{

///simple check to see if file exists
bool FileExists(const char *fname);

//-- Read routines

/// \name Reading routines
/// see \ref io.cxx for implementation
//@{
///Reads the header information
Int_t ReadHeader(Options &opt);
///Reads halo particle id data
HaloTreeData *ReadData(Options &opt);
///Reads individual halo data file
HaloData *ReadHaloData(string &fname, Int_t &nhalos);
///Reads NIFTY AHF individual halo data file
HaloData *ReadNIFTYData(string &fname, Int_t &nhalos, int idcorrectflag=0, int hidoffset=1);
///Reads a Void catalog file
HaloData *ReadVoidData(string &fname, Int_t &nhalos, int idcorrectflag=0, int hidoffset=1);
///Reads VELOCIraptor like Group Catalog Data. Can adjust so that only particles of some type are check for cross matching
HaloData *ReadHaloGroupCatalogData(string &infile, Int_t &numhalos, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH, int iverbose=0);
//@}

/// \name routines used to read VELOCIraptor output
/// see \ref stfio.cxx for implementation
//@{
///minor interface routine associated with reading VELOCIraptor output, openfiles appropriate files
void OpenBinaryorAsciiFiles(string &infile, int ibinary, int numfiletypes, int k, int mpi_ninput, int ifieldhalos, int itypematch, fstream &Fgroup, fstream &Fpart, fstream &Fupart, fstream &Fsgroup, fstream &Fspart, fstream &Fsupart, fstream &Fparttype, fstream &Fuparttype, fstream &Fsparttype, fstream &Fsuparttype, int iverbose=0);
///minor interface for closing appropriate files
void CloseBinaryorAsciiFiles(fstream &Fgroup, fstream &Fpart, fstream &Fupart, fstream &Fsgroup, fstream &Fspart, fstream &Fsupart, fstream &Fparttype, fstream &Fuparttype, fstream &Fsparttype, fstream &Fsuparttype, int ifieldhalos, int itypematch);
#ifdef USEHDF
void OpenHDFFiles(string &infile, int numfiletypes, int k, int mpi_ninput, int ifieldhalos, int itypematch, H5File &Fgroup, H5File &Fpart, H5File &Fupart, H5File &Fsgroup, H5File &Fspart, H5File &Fsupart, H5File &Fparttype, H5File &Fuparttype, H5File &Fsparttype, H5File &Fsuparttype, int iverbose=0);
void CloseHDFFiles(H5File &Fgroup, H5File &Fpart, H5File &Fupart, H5File &Fsgroup, H5File &Fspart, H5File &Fsupart, H5File &Fparttype, H5File &Fuparttype, H5File &Fsparttype, H5File &Fsuparttype, int ifieldhalos, int itypematch);
#endif
///read routine to load number of files the VELOCIraptor output is split between
inline void STFReadNumFileInfoAndCorrectNumFile(int &itask, int &nprocs, int &nmpicount, int &mpi_ninput, fstream &Fgroup, fstream &Fsgroup,
#ifdef USEHDF
    H5File &Fhdfgroup, DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary, int ifieldhalos);
///read information from the group catalog file and correct the number of files if necessary
inline void STFReadNumGroups(unsigned long &nglocal, unsigned long &TotalNumberofHalos, unsigned long &nsglocal, fstream &Fgroup, fstream &Fsgroup,
#ifdef USEHDF
    H5File &Fhdfgroup, H5File &Fhdfsgroup, DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary, int ifieldhalos);
inline void STFReadNumData(unsigned long &nids, unsigned long &nsids, unsigned long &nuids, unsigned long &nsuids, 
    unsigned long &nidstot, unsigned long &nsidstot, unsigned long &nuidstot, unsigned long &nsuidstot, 
    fstream &Fpart, fstream &Fupart, fstream &Fspart, fstream &Fsupart, 
    fstream &Fparttype, fstream &Fuparttype, fstream &Fsparttype, fstream &Fsuparttype, 
#ifdef USEHDF
    H5File &Fhdfpart, H5File &Fhdfspart, H5File &Fhdfupart, H5File &Fhdfsupart, 
    H5File &Fhdfparttype, H5File &Fhdfsparttype, H5File &Fhdfuparttype, H5File &Fhdfsuparttype, 
    DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary, int ifieldhalos, int itypematch);
//@}


/// \name Write routines
/// see \ref io.cxx for implementation
//@{
///Writes the merger tree
void WriteHaloMergerTree(Options &opt, ProgenitorData **p, HaloTreeData *h);
///writes a full graph
void WriteHaloGraph(Options &opt, ProgenitorData **p, DescendantData **d, HaloTreeData *h);
///writes a cross comparison between two catalogs
void WriteCrossComp(Options &opt, ProgenitorData **p, HaloTreeData *h);
//@}

#ifdef USEMPI
/// \name  for mpi related reads
/// see \ref io.cxx, \ref stfio.cxx, \ref mpiroutines.cxx for implementation
//@{
///Reads the number of haloes in the files, wrapper routine
Int_t ReadNumberofHalos(Options &opt, Int_t *numhalos);
///Reads the number of particle in haloes in the files, wrapper routine
Int_t ReadNumberofParticlesInHalos(Options &opt, Int_t *numpartinhalos);
///Reads Number of halos from VELOCIraptor Group catalogs
Int_t MPIReadHaloGroupCatalogDataNum(string &infile, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH);
///Reads Number of particles in halos from VELOCIraptor Group catalogs
Int_t MPIReadHaloGroupCatalogDataParticleNum(string &infile, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH);
///Reads data to allocate memory, useful for mpi
HaloData *MPIReadHaloGroupCatalogDataAllocation(string &infile, Int_t &numhalos, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH);
///Reads VELOCIraptor like Group Catalog Data with memory already allocated for MPI version. 
void MPIReadHaloGroupCatalogData(string &infile, Int_t &numhalos, HaloData *&Halo, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH, int iverbose=0);
///load balance input data
void MPILoadBalanceSnapshots(Options &);
///as load balancing is long, read/write the load balance
void MPIWriteLoadBalance(Options &);
///as load balancing is long, read/write the load balance
int MPIReadLoadBalance(Options &);
//@}

///\name  for mpi related meshing of data
/// see \ref mpiroutines.cxx for implementation
//@{
void MPIUpdateProgenitorsUsingDescendants(Options &opt, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen);
//@}
#endif

//@}

/// \name Subroutines calculates merits between 
/// see \ref crossmatch.cxx for implementation
//@{

/// determine the cross matches of halos in h1 in "progenitor list" h2
/// the routine also stores whether the progenitor list has been updated through the ilistupdated variable.
/// the routine can also be provided a reference list and the time step of the current halo list used to identify progenitors
/// and will only search current halo list if the reference list doesn't meet the criteria for viable progenitors
/// which is currently whether there are any in the reference list
ProgenitorData *CrossMatch(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData* &h2, unsigned int*&pfof2, int &ilistupdated, int istepval=1, ProgenitorData *refprogen=NULL);
///clean cross matches of duplicate entries so that a progenitor can have ONLY ONE descendent. 
void CleanCrossMatch(const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, ProgenitorData *&pprogen);
///fill in empty links of the reference list with another progenitor list produced using same reference snapshot but different linking snapshot. 
///Allows for multiple steps in snapshots to be used. 
void UpdateRefProgenitors(Options &opt, const Int_t numhalos,ProgenitorData *&pref, ProgenitorData *&ptemp, DescendantDataProgenBased **&pprogendescen, Int_t itimedescen);
///similar to \ref CrossMatch but for descendants
DescendantData *CrossMatchDescendant(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData* &h2, unsigned int*&pfof2, int &ilistupdated, int istepval=1, DescendantData *refdescen=NULL);
///similar to \ref CleanCrossMatch but for descendants
void CleanCrossMatchDescendant(const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, DescendantData *&pdescen);
///similar to \ref UpdateRefProgenitors but for descendants
void UpdateRefDescendants(Options &opt, const Int_t numhalos,DescendantData *&pref, DescendantData *&ptemp);
///builds the possible set of descendents using candidate progenitors
void BuildProgenitorBasedDescendantList(Int_t itimeprogen, Int_t itimedescen, Int_t nhalos, ProgenitorData *&pprogen, DescendantDataProgenBased **&pprogendescen, int istep=1);
/// removes descendant links of individual halo 
void RemoveLinksProgenitorBasedDescendantList(Int_t itimedescen, Int_t ihaloindex, ProgenitorData &pprogen, DescendantDataProgenBased **&pprogendescen);
///Cleans up progenitor list using candidate descendent list build using progenitor search. Ensures objects ony have a single descendant
void CleanProgenitorsUsingDescendants(Int_t i, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen);

//@}

/// \name for mapping ids to index routines
/// see \ref crossmatch.cxx for implementation
//@{
///map particle id to index position
void MapPIDStoIndex(Options &opt, HaloTreeData *&pht);
///make sure particle ids are acceptable values for generating links
void IDcheck(Options &opt,HaloTreeData *&pht);

void simplemap(long unsigned &i);

//adjust the mappable halo ids by adding a temporally unique value to each id
//useful to ensure that ids between VELOCIraptor and the tree match
void UpdateHaloIDs(Options &opt, HaloTreeData *&pht);

//@}

/// \name check particle types when building list that will be cross matched
//@{
///map particle id to index position
int CheckType(unsigned int t, int tmatch);
//@}

/// \name Extra routines
/// see \ref crossmatch.cxx for implementation
//@{

///build group size array
Int_t *BuildNumInGroup(const Int_t nbodies, const Int_t numgroups, Int_t *pfof);
///build array such that array is pglist[group][]={particle list}
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof);
///build the group particle index list (doesn't assume particles are in ID order and stores index of particle)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *Part);
///build pglist but doesn't assume particles are in ID order
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Int_t *ids);
///build the group particle arrays need for unbinding procedure
Particle **BuildPartList(const Int_t numgroups, Int_t *numingroup, Int_t **pglist, Particle* Part);
///build a particle list subset using array of indices
Particle *BuildPart(Int_t numingroup, Int_t *pglist, Particle* Part);
///build the Head array which points to the head of the group a particle belongs to
Int_t *BuildHeadArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///build the Next array which points to the next particle in the group
Int_t *BuildNextArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///build the Len array which stores the length of the group a particle belongs to 
Int_t *BuildLenArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);
///build the GroupTail array which stores the Tail of a group
Int_t *BuildGroupTailArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist);

//@}

/// \name Extra utility routines
/// see \ref utilities.cxx for implementation
//@{

///Get current time in Milliseconds
int GetMilliCount();
///Get span in milliseconds
int GetMillSpan(int );
///get a time
double MyGetTime();
//@}

#endif


