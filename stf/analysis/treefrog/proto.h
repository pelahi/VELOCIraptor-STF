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
void GetParamFile(Options &opt);
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
///writes descendant halo merger tree
void WriteHaloMergerTree(Options &opt, DescendantData **p, HaloTreeData *h);
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
unsigned long ReadNumberofHalos(Options &opt, unsigned long *numhalos);
///Reads the number of particle in haloes in the files, wrapper routine
unsigned long ReadNumberofParticlesInHalos(Options &opt, unsigned long *numpartinhalos);
///Reads Number of halos from VELOCIraptor Group catalogs
unsigned long MPIReadHaloGroupCatalogDataNum(string &infile, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH);
///Reads Number of particles in halos from VELOCIraptor Group catalogs
unsigned long MPIReadHaloGroupCatalogDataParticleNum(string &infile, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH);
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
///update descendant based progenitor information across mpi domain
void MPIUpdateProgenitorsUsingDescendants(Options &opt, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen);
///update progenitor based descendant information across mpi domain
void MPIUpdateDescendantUsingProgenitors(Options &opt, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen);
///update descendant information across mpi domain
void MPIUpdateDescendants(Options &opt, HaloTreeData *&pht, DescendantData **&pdescen);

///helper routines to send information
void MPISendProgenitorsUsingDescendants(int recvtask, int isnap, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen);
///helper routines to recv and process information
void MPIRecvProgenitorsUsingDescendants(int sendtask, int isnap, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen);

///helper routines to send information
void MPISendDescendantsUsingProgenitors(int recvtask, int isnap,HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen);
///helper routines to recv and process information
void MPIRecvDescendantsUsingProgenitors(int sendtask, int isnap,HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen);


///helper routines to send information
void MPISendDescendants(int recvtask, int isnap,HaloTreeData *&pht, DescendantData **&pdescen);
///helper routines to recv and process information
void MPIRecvDescendants(int sendtask, int isnap,HaloTreeData *&pht, DescendantData **&pdescen);
//@}
#endif

//@}

/// \name Subroutines calculates merits between
/// see \ref crosscheck.cxx for implementation
//@{

///Calculate merit between a match
Double_t CalculateMerit(Options &opt, UInt_t n1, UInt_t n2, HaloData &h1, HaloData &h2, UInt_t hindex=0,UInt_t *sharepartlist=NULL, Double_t *ranking2=NULL);

/// determine the cross matches of halos in h1 in "progenitor list" h2
/// the routine also stores whether the progenitor list has been updated through the ilistupdated variable.
/// the routine can also be provided a reference list and the time step of the current halo list used to identify progenitors
/// and will only search current halo list if the reference list doesn't meet the criteria for viable progenitors
/// which is currently whether there are any in the reference list
ProgenitorData *CrossMatch(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData* &h2, unsigned int*&pfof2, int &ilistupdated, int istepval=1, ProgenitorData *refprogen=NULL);
///get Progenitor match for individual object, return if match found
int CrossMatchProgenitorIndividual(Options &opt, Int_t i,
    const long unsigned nhalos1, const long unsigned nhalos2,
    HaloData *&h1, HaloData *&h2,
    unsigned int *&pfof2,
    int istepval,
    ProgenitorData *&p1,
    unsigned int *&sharelist,
    unsigned int *&halolist,
    long unsigned offset
    //unsigned int *&sharepartlist,
    //unsigned int *&pranking2,
    //Double_t *&rankingsum
);

///clean cross matches of duplicate entries so that a progenitor can have ONLY ONE descendent.
void CleanCrossMatch(const int istepval, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, ProgenitorData *&pprogen);
///fill in empty links of the reference list with another progenitor list produced using same reference snapshot but different linking snapshot.
///Allows for multiple steps in snapshots to be used.
void UpdateRefProgenitors(Options &opt, const Int_t numhalos,ProgenitorData *&pref, ProgenitorData *&ptemp, DescendantDataProgenBased **&pprogendescen, Int_t itime);
///builds the possible set of descendents using candidate progenitors
void BuildProgenitorBasedDescendantList(Int_t itimeprogen, Int_t itimedescen, Int_t nhalos, ProgenitorData *&pprogen, DescendantDataProgenBased **&pprogendescen, int istep=1);
/// removes descendant links of individual halo
void RemoveLinksProgenitorBasedDescendantList(Int_t itimedescen, Int_t ihaloindex, ProgenitorData &pprogen, DescendantDataProgenBased **&pprogendescen);
///Cleans up progenitor list using candidate descendent list build using progenitor search. Ensures objects ony have a single descendant
void CleanProgenitorsUsingDescendants(Int_t i, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen, int iopttemporalmerittype);

///similar to \ref CrossMatch but for descendants
DescendantData *CrossMatchDescendant(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData* &h2, unsigned int*&pfof2, int &ilistupdated, int istepval=1, unsigned int *ranking2=0, DescendantData *refdescen=NULL);
///get descendant match for individual object, return if match found
int CrossMatchDescendantIndividual(Options &opt, Int_t i,
    const long unsigned nhalos1, const long unsigned nhalos2,
    HaloData *&h1, HaloData *&h2,
    unsigned int *&pfof2,
    int istepval, int initdtopval,
    DescendantData *&d1,
    unsigned int *&sharelist,
    unsigned int *&halolist,
    long unsigned offset, long unsigned offset2,
    unsigned int *&sharepartlist,
    unsigned int *&pranking2,
    Double_t *&rankingsum
);

///updates the haloids stored in the descendant list
void UpdateDescendantIndexing(const int istepval, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, DescendantData *&p1);
///prunes the descendant list
void CleanCrossMatchDescendant(Options &opt, Int_t itime, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdescen);

///builds the possible set of progenitors using candidate descendants, similar to \ref BuildProgenitorBasedDescendantList.
void BuildDescendantBasedProgenitorList(Int_t itimedescen, Int_t nhalos, DescendantData *&pdecen, ProgenitorDataDescenBased *&pdescenprogen, int istep=1);
///updates descendants data structure so that an object knows what type of progenitor it is to its descendants. Specifically this finds all progenitors of
///a halo looking back a time istep and ranks their descendant to progenitor match based on their merit (a temporally local ranking)
void UpdateDescendantUsingDescendantBasedProgenitorList(Int_t nhalos, DescendantData *&pdescen, ProgenitorDataDescenBased *&pdescenprogen, int istep, Double_t meritlimit);
/// Updates the descandant data structure with new information and removes/updates old links
void UpdateRefDescendants(Options &opt, const Int_t numhalos,DescendantData *&pref, DescendantData *&ptemp, ProgenitorDataDescenBased **&pdescenprogen, Int_t itime);
///removes progenitor links of individual halo stored in the \ref ProgenitorDataDescenBased structure.
void RemoveLinksDescendantBasedProgenitorList(Int_t itime, Int_t ihaloindex, DescendantData &pdescen, ProgenitorDataDescenBased **&pdescenprogen);
///Adds progenitor links of individual halo stored in the \ref ProgenitorDataDescenBased structure.
void AddLinksDescendantBasedProgenitorList(Int_t itime, Int_t ihaloindex, DescendantData &pdescen, ProgenitorDataDescenBased **&pdescenprogen);

///Rank the progenitors in the descendant tree using the descendants and a general temporal merit.
void RankDescendantProgenitors(Int_t i, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdecen, int iopttemporalmerittype);
///clean the descendant list to adjust descedant rankings to minimize number of objects with no primary progenitors
void CleanDescendantsForMissingProgenitors(Options &opt, Int_t itime, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdescen);

///Reranks descendants based on descendant to progenitor ranking and then merit.
void RerankDescendants(Options &opt, HaloTreeData *&pht, DescendantData **&pdescen);


//@}

/// \name for mapping ids to index routines
/// see \ref idroutines.cxx for implementation
//@{
///wrapper routine that adjusts data to use memory efficient maps
void MemoryEfficientMap(Options &opt,HaloTreeData *&pht);
///generate a map for particle ids to index and store it in the set, which can be searched to
///relate pid to index
map<IDTYPE, IDTYPE> ConstructMemoryEfficientPIDStoIndexMap(Options &opt, HaloTreeData *&pht);
///save the particle id to index map
void SavePIDStoIndexMap(Options &,map<IDTYPE, IDTYPE>&);
///read the particle id to index map from file
int ReadPIDStoIndexMap(Options &, map<IDTYPE, IDTYPE>&);

///map particle id to index position
void MapPIDStoIndex(Options &opt, HaloTreeData *&pht, map<IDTYPE, IDTYPE> &);
///map particle id to index position
void MapPIDStoIndex(Options &opt, HaloTreeData *&pht);

///make sure particle ids are acceptable values for generating links
void IDcheck(Options &opt,HaloTreeData *&pht);
///simple mapping function
void simplemap(IDTYPE &i);

///adjust the mappable halo ids by adding a temporally unique value to each id
///useful to ensure that ids between VELOCIraptor and the tree match
void UpdateHaloIDs(Options &opt, HaloTreeData *&pht);

//@}

/// \name check particle types when building list that will be cross matched
//@{
///map particle id to index position
int CheckType(unsigned int t, int tmatch);
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
