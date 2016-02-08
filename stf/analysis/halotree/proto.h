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

///Reads the header information
Int_t ReadHeader(Options &opt);
///Reads halo particle id data
HaloTreeData *ReadData(Options &opt);
///Reads individual halo data file
HaloData *ReadHaloData(char *fname, Int_t &nhalos);
///Reads NIFTY AHF individual halo data file
HaloData *ReadNIFTYData(char *fname, Int_t &nhalos, int idcorrectflag=0);
///Reads VELOCIraptor like Group Catalog Data. Can adjust so that only particles of some type are check for cross matching
HaloData *ReadHaloGroupCatalogData(char* infile, Int_t &numhalos, int mpi_ninput=0, int ibinary=1, int ifieldhalos=1, int itypesort=ALLTYPEMATCH);

///Writes the merger tree
void WriteHaloMergerTree(Options &opt, ProgenitorData **p, HaloTreeData *h);
void WriteHaloGraph(Options &opt, ProgenitorData **p, DescendantData **d, HaloTreeData *h);
///writes a cross comparison between two catalogs
void WriteCrossComp(Options &opt, ProgenitorData **p, HaloTreeData *h);

//@}

/// \name Subroutines calculates merits between 
/// see \ref crossmatch.cxx for implementation
//@{

/// determine the cross matches of halos in h1 in "progenitor list" h2
ProgenitorData *CrossMatch(Options &opt, const Int_t nbodies, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData* &h2, long unsigned *&pglist, long unsigned *&noffset, long unsigned *&pfof2, int istepval=1);
///clean cross matches of duplicate entries so that a progenitor can have ONLY ONE descendent. 
void CleanCrossMatch(const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, ProgenitorData *&pprogen);
///fill in empty links of the reference list with another progenitor list produced using same reference snapshot but different linking snapshot. 
///Allows for multiple steps in snapshots to be used. 
void UpdateRefProgenitors(const Int_t numhalos,ProgenitorData *&pref, ProgenitorData *&ptemp);
///similar to \ref CrossMatch but for descendants
DescendantData *CrossMatchDescendant(Options &opt, const Int_t nbodies, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData* &h2, long unsigned *&pglist, long unsigned *&noffset, long unsigned *&pfof2, int istepval=1);
///similar to \ref CleanCrossMatch but for descendants
void CleanCrossMatchDescendant(const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, DescendantData *&pdescen);
///similar to \ref UpdateRefProgenitors but for descendants
void UpdateRefDescendants(const Int_t numhalos,DescendantData *&pref, DescendantData *&ptemp);

//@}

/// \name for mapping ids to index routines
//@{
///map particle id to index position
void MapPIDStoIndex(Options &opt, HaloTreeData *&pht);

void simplemap(long unsigned &i);

//@}

/// \name check particle types when building list that will be cross matched
//@{
///map particle id to index position
int CheckType(unsigned int t, int tmatch);
//@}

/// \name Extra routines
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

#endif


