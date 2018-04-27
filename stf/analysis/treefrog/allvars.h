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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/timeb.h>

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

#ifdef USEOPENMP
#include <omp.h>
#endif

#ifdef USEMPI
#include <mpi.h>
///Includes external global variables used for MPI version of code
#include "mpivar.h"
#endif

///if using HDF API
#ifdef USEHDF
#include "H5Cpp.h"
using namespace H5;
#endif

///if using ADIOS API
#ifdef USEADIOS
#include "adios.h"
#endif

using namespace std;
using namespace Math;
using namespace NBody;

//-- Structures and external variables


/// \name Specific structure type, allow for other types beside HALO
/// \todo will have to alter this to allow for subvois possibly even subfilaments (though I'm not certain what the hell this would be)
//@{
#define HALOSTYPE 10
#define WALLSTYPE 2
#define VOIDSTYPE 1
#define FILAMENTSTYPE 3
//@}

/// \name Cross-match types
//@{

//subsets made
///call \ref NsharedN1N2
#define MERITNsharedN1N2 1
#define MERITNsharedN1 2
#define MERITNshared 3
#define MERITNsharedcombo 4
#define MERITRankWeighted 5
#define MERITRankWeightedBoth 6
//@}

/// \name Input format
//@{
#define DSUSSING 1
#define DCATALOG 2
#define DNIFTY 3
#define DVOID 4
//@}

/// \name operation/type of catalog
//@{
#define DTREE 0
#define DCROSSCAT 1
#define DGRAPH 2
#define DDESCENDANT 3
//@}

/// \name defining direction of search when constructing a tree
//@{
///standard progenitor search
#define SEARCHPROGEN 0
///standard descendant search
#define SEARCHDESCEN 1
///search both ways
#define SEARCHALL -1
//@}

/// \name defining the type of core particle list matching
//@{
///standard search
#define PARTLISTNOCORE 0
///limit search to core
#define PARTLISTCORE 1
///limit search to core in both searched particles and matched particles to identify best match after finding all matches
#define PARTLISTCORECORE 2
///limit search to core in both searched particles and matched particles ONLY
#define PARTLISTCORECOREONLY 3
//@}

/// \name Mapping functions
//@{
#define DNOMAP 0
#define DSIMPLEMAP 1
#define DMEMEFFICIENTMAP -1
//@}


/// \name nifty format particle type for subselecting particle matching
#define NIFTYDMTYPE 1

/// \name defining particle types for subselecting particle matching when loading velociraptor group catalogs
//@{
#define ALLTYPEMATCH -1
#define GASTYPEMATCH 0
#define DMTYPEMATCH 1
#define DMGASTYPEMATCH -2
#define STARTYPEMATCH 4
//@}

/// \name input formats for velociraptor output
//@{
#define INBINARY 1
#define INHDF 2
#define INADIOS 3
#define INASCII 0
//@}

/// \name output formats
//@{
#define OUTASCII 0
#define OUTBINARY 1
#define OUTHDF 2
//@}

/// \name defining types of multisnapshot linking done
//@{
///missing link
#define MSLCMISSING 0
///higher merit
#define MSLCMERIT 1
///higher primary progenitor
#define MSLCPRIMARYPROGEN 2
///higher merit and primary progenitor
#define MSLCMERITPRIMARYPROGEN 3
//@}

/// \name defining types of optimal temporal merit criteria
//@{
/// simple temporal merit
#define GENERALIZEDMERITTIME 0
/// temporal merit and also as close to the primary progenitor as possible
#define GENERALIZEDMERITTIMEPROGEN 1
//@}

/// \name parameters for the temporal merit function
//@{
///the index of the temporal weight in delta t ^(-alpha)
#define ALPHADELTAT 0.12
//@}

///VALUE USED TO MAKE UNIQUE halo identifiers
#ifdef LONGINT
#define HALOIDSNVAL 1000000000000L
#else
#define HALOIDSNVAL 1000000
#endif

/// \name OpenMP parameters for load balancing
#ifdef USEOPENMP
#define OMPCHUNKSIZE 100000UL
#endif

/// \name MPI parameters for load balancing
#ifdef USEMPI
///particle based splitting
#define MPIPARTICLEBALANCE 0
///halo based splitting
#define MPIHALOBALANCE 1
#endif

#ifdef TREEFROGLONGIDS
typedef long long IDTYPE;
#elif defined(TREEFROGLONGUIDS)
typedef long unsigned IDTYPE;
#elif defined(TREEFROGINTIDS)
typedef int IDTYPE;
#else
typedef int unsigned IDTYPE;
#endif

/// Options structure stores useful variables that have user determined values which are altered by \ref GetArgs in \ref ui.cxx
struct Options
{
    ///\name filenames
    //@{
    char *fname,*outname, *configname;
    //@}
    ///number of snapshots
    int numsnapshots;
    ///nubmer of steps integrated over to find links (two steps needs three snapshots)
    int numsteps;
    ///maximum id value, used to allocate an array of this size so that ids can be mapped to an index and thus easily accessible.
    long unsigned MaxIDValue;
    ///total number of haloes across all snapshots
    long unsigned TotalNumberofHalos;
    ///read halo positions
    int iposinfoflag;

    /// store description of code
    string description;

    ///type of search, going back in time (progenitor searching) or going forward in time (descendant searching) or both
    int isearchdirection;
    ///type of merit
    int imerittype;
    ///cross match shared particle number significance, that is match only when quantity is above mlsig*some measure of noise, here defined as
    ///\$f N_2^{1/2} \$f
    Double_t mlsig;
    ///cross match merit limit for deciding whether to search previous snapshot
    Double_t meritlimit;
    ///allowed merit ratio when comparing merits between two viable connections for swapping
    Double_t meritratiolimit;
    ///merit ratio when deciding that have significant ambiguity no ideal descendant/progenitor found
    Double_t meritratioambiguitylimit;

    //if particle ids are exclusive to a structure
    int iexclusiveids;

    ///particle type cross matching subselection
    int itypematch;
    ///when using multiple links, how links should be updated, either only for those missing links, or better merit found
    int imultsteplinkcrit;
    ///set the optimal temporal merit
    int iopttemporalmerittype;

    ///flag for whether default merger tree produced, a simple cross comparsion between two catalogs is run or a full graph is constructed
    int icatalog;

    ///output format
    int outputformat;
    ///format of data written
    int outdataformat;

    ///\name io format flags for basic format, whether data split into multiple files, binary, and field objects separate
    //@{
    int ioformat,nmpifiles,ibinary,ifield;
    //@}

    //for mapping particle ids to index
    //@{
    int imapping;
    void(*mappingfunc)(IDTYPE &);
    //@}

    ///for offseting the snapshot id value when writing to a file
    Int_t snapshotvaloffset;

    ///for adjusting halo ids by (current_snap + \ref snapshotvaloffset )* times this value
    long long haloidval;

    ///for adjusting halo ids by a simple offset
    Int_t haloidoffset;

    ///whether to use core particles only or not, 0 is use all particles to identify matches, PARTLISTCORE is only core particles defined by particle fraction, and PARTLISTCORECORE only match
    ///core particles to other core particles
    int icorematchtype;
    ///to adjust the number of particles used to define a merit. Note that to be of use, particles should be in some type of order
    ///for instance binding energy order. VELOCIraptor outputs particles in decreasing binding energy order
    ///so using the first particle_frac of a halo means using the most bound set of particles
    Double_t particle_frac;
    ///minimum number of particles used to derive merit.
    Int_t min_numpart;
    ///maximum number of particles used to derive merit. If -1, no maximum
    Int_t max_numpart;

    ///to fix ids for nifty project
    int idcorrectflag;

    ///verbose output flag
    int iverbose;
    ///use default optimized values, ignore some input values
    int idefaultvalues;

#ifdef USEMPI
    ///number of items (halos or particles in halos) across various snapshots desired. Used for load balancing
    Int_t numpermpi;
    ///number of expected mpi threads, useful for running code with a single thread to determine load balancing
    ///which then writes the mpi load balancing file and exists
    int ndesiredmpithreads;
    ///whether mpi threads write in parallel
    int iwriteparallel;
    ///Set if the load balancing is either halo or particle based splitting
    int impiloadbalancesplitting;

#endif

    Options()
    {
        fname=outname=NULL;

        numsnapshots=2;
        numsteps=1;
        MaxIDValue=512*512*512;
        TotalNumberofHalos=0;
        iexclusiveids=1;
        iposinfoflag=0;

        isearchdirection=SEARCHPROGEN;
        imerittype=MERITNsharedN1N2;
        mlsig=0.1;
        meritlimit=0.05;
        itypematch=ALLTYPEMATCH;
        imultsteplinkcrit=MSLCMERIT;
        iopttemporalmerittype=GENERALIZEDMERITTIMEPROGEN;

        ioformat=DCATALOG;
        icatalog=DTREE;

        nmpifiles=0;
        ibinary=1;
        ifield=1;

        imapping=DNOMAP;//no mapping
        mappingfunc=NULL;

        snapshotvaloffset=0;
        haloidval=0;
        idcorrectflag=0;
        outputformat=OUTASCII;
        outdataformat=0;
        haloidoffset=0;

        icorematchtype=PARTLISTNOCORE;
        particle_frac=0.2;
        min_numpart=20;
        max_numpart=-1;
        meritratiolimit=4.0;
        meritratioambiguitylimit=1.5;

        iverbose=1;
        idefaultvalues=1;
#ifdef USEMPI
        numpermpi=0;
        ndesiredmpithreads=0;
        iwriteparallel=0;
        impiloadbalancesplitting=MPIHALOBALANCE;
#endif
    }
};

/*!
    Halo data structure
*/
struct HaloData{
    long unsigned haloID;
    long unsigned NumberofParticles;
    IDTYPE *ParticleID;
    float Xcm[3], Vcm[3];
    float Rmax, Vmax;
    HaloData(long unsigned np=0){
        NumberofParticles=np;
        ParticleID=NULL;
        if (NumberofParticles>0) {
            ParticleID=new IDTYPE[NumberofParticles];
        }
    }
    void Alloc(long unsigned np=0){
        //if (NumberofParticles>0){
        if (ParticleID!=NULL){
            delete[] ParticleID;
            ParticleID=NULL;
        }
        NumberofParticles=np;
        if (NumberofParticles>0) {
            ParticleID=new IDTYPE[NumberofParticles];
        }
    }
    ~HaloData(){
        //if (NumberofParticles>0){
        if (ParticleID!=NULL){
            delete[] ParticleID;
        }
    }
};

struct HaloTreeData{
    Int_t numhalos;
    HaloData *Halo;
    HaloTreeData(){
        numhalos=0;
        Halo=NULL;
    }
    ~HaloTreeData(){
        if (numhalos>0) delete[] Halo;
    }
};

/*!
    Structure used to keep track of a structure's progenitor/parent structure using temporal/spatial information
*/
struct ProgenitorData
{
    ///structure type and number of links (progenitors)
    //@{
    int stype;
    int NumberofProgenitors;
    //@}
    ///store list of progenitors
    long unsigned* ProgenitorList;
    ///store the merit value
    float *Merit;
    ///store the fraction of shared particles
    float *nsharedfrac;
    ///store number of steps back in time progenitor found
    int istep;

    ProgenitorData(){
        NumberofProgenitors=0;
        ProgenitorList=NULL;
        Merit=NULL;
        nsharedfrac=NULL;
        istep=0;
    }
    ~ProgenitorData(){
        if (NumberofProgenitors>0) {
            delete[] ProgenitorList;
            delete[] Merit;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
    }
    ProgenitorData &operator=(const ProgenitorData &p){
        if (NumberofProgenitors>0) {
            delete[] ProgenitorList;
            delete[] Merit;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
        NumberofProgenitors=p.NumberofProgenitors;
        ProgenitorList=new long unsigned[NumberofProgenitors];
        Merit=new float[NumberofProgenitors];
        for (int i=0;i<NumberofProgenitors;i++) {
            ProgenitorList[i]=p.ProgenitorList[i];
            Merit[i]=p.Merit[i];
        }
        if (p.nsharedfrac!=NULL) {
            nsharedfrac=new float[NumberofProgenitors];
            for (int i=0;i<NumberofProgenitors;i++) nsharedfrac[i]=p.nsharedfrac[i];
        }
        istep=p.istep;
    }

};

/*!
    Structure used to keep track of a structure's descendant/child structure using temporal/spatial information
*/
struct DescendantData
{
    ///structure type and number of links (descendants)
    //@{
    int stype;
    int NumberofDescendants;
    //@}
    ///store list of descendants
    long unsigned* DescendantList;
    ///store the merit value
    float *Merit;
    ///store the fraction of shared particles
    float *nsharedfrac;
    ///store type of descendant progenitor relation, that is whether object is the primary progenitor of its descendant or not
    unsigned short *dtoptype;

    ///store number of steps forward in time descedant found
    int istep;

    DescendantData(){
        NumberofDescendants=0;
        DescendantList=NULL;
        Merit=NULL;
        nsharedfrac=NULL;
        dtoptype=NULL;
        istep=0;
    }
    DescendantData(int &nd, long unsigned *dl, float *m, unsigned short *dtop, int &step, float *ns=NULL){
        NumberofDescendants=nd;
        DescendantList=new long unsigned[NumberofDescendants];
        Merit=new float[NumberofDescendants];
        for (int i=0;i<NumberofDescendants;i++) {
            DescendantList[i]=dl[i];
            Merit[i]=m[i];
        }
        dtoptype=new unsigned short[NumberofDescendants];
        for (int i=0;i<NumberofDescendants;i++) dtoptype[i]=dtop[i];
        istep=step;
        if (ns!=NULL) {
            nsharedfrac=new float[NumberofDescendants];
            for (int i=0;i<NumberofDescendants;i++) nsharedfrac[i]=ns[i];
        }
    }
    ~DescendantData(){
        if (NumberofDescendants>0) {
            delete[] DescendantList;
            delete[] Merit;
            if (dtoptype!=NULL) delete[] dtoptype;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
    }
    void Set(int &nd, long unsigned *dl, float *m, unsigned short *dtop, int &step, float *ns=NULL){
        if (NumberofDescendants>0) {
            delete[] DescendantList;
            delete[] Merit;
            if (dtoptype!=NULL) delete[] dtoptype;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
        NumberofDescendants=nd;
        DescendantList=new long unsigned[NumberofDescendants];
        Merit=new float[NumberofDescendants];
        for (int i=0;i<NumberofDescendants;i++) {
            DescendantList[i]=dl[i];
            Merit[i]=m[i];
        }
        dtoptype=new unsigned short[NumberofDescendants];
        for (int i=0;i<NumberofDescendants;i++) dtoptype[i]=dtop[i];
        istep=step;
        if (ns!=NULL) {
            nsharedfrac=new float[NumberofDescendants];
            for (int i=0;i<NumberofDescendants;i++) nsharedfrac[i]=ns[i];
        }
    }

    DescendantData &operator=(const DescendantData &d){
        if (NumberofDescendants>0) {
            delete[] DescendantList;
            delete[] Merit;
            if (dtoptype!=NULL) delete[] dtoptype;
            if (nsharedfrac!=NULL) delete[] nsharedfrac;
        }
        NumberofDescendants=d.NumberofDescendants;
        DescendantList=new long unsigned[NumberofDescendants];
        Merit=new float[NumberofDescendants];
        for (int i=0;i<NumberofDescendants;i++) {
            DescendantList[i]=d.DescendantList[i];
            Merit[i]=d.Merit[i];
        }
        if (d.dtoptype!=NULL) {
            dtoptype=new unsigned short[NumberofDescendants];
            for (int i=0;i<NumberofDescendants;i++) dtoptype[i]=d.dtoptype[i];
        }
        if (d.nsharedfrac!=NULL) {
            nsharedfrac=new float[NumberofDescendants];
            for (int i=0;i<NumberofDescendants;i++) nsharedfrac[i]=d.nsharedfrac[i];
        }
        istep=d.istep;
    }

};

/*!
    Structure used to keep track of a structure's candidate descendants based on the progenitor construction
*/
struct DescendantDataProgenBased
{
    ///structure type and number of links (descendants)
    //@{
    int stype;
    int NumberofDescendants;
    //@}

    ///store list of descendants in the form of halo index and temporal index
    vector<long unsigned> haloindex;
    vector<int unsigned> halotemporalindex;
    ///store the merit value
    vector<float> Merit;
    //store the integer time diff to descendant
    vector<int> deltat;
    //store descendant type, primary or not primary
    vector<int> descentype;
#ifdef USEMPI
    //store local number of mpi items
    int nlocal;
    //store which task this halo progenitor is located on
    vector<int> MPITask;
    //if using mpi, need to keep track of items removed that
    //stuff can be properly mergered
    vector<long unsigned> removalhaloindex;
    vector<int unsigned> removalhalotemporalindex;
#endif
    DescendantDataProgenBased(int reservesize=4){
        NumberofDescendants=0;
        //reserve some space so mimimize number of reallocations
        haloindex.reserve(reservesize);
        halotemporalindex.reserve(reservesize);
        Merit.reserve(reservesize);
        deltat.reserve(reservesize);
        descentype.reserve(reservesize);
#ifdef USEMPI
        MPITask.reserve(reservesize);
        removalhaloindex.reserve(reservesize);
        removalhalotemporalindex.reserve(reservesize);
#endif
    }
    ~DescendantDataProgenBased(){
    }
    DescendantDataProgenBased &operator=(const DescendantDataProgenBased &d){
        NumberofDescendants=d.NumberofDescendants;
        haloindex.resize(NumberofDescendants);
        halotemporalindex.resize(NumberofDescendants);
        Merit.resize(NumberofDescendants);
        descentype.reserve(NumberofDescendants);
        deltat.resize(NumberofDescendants);
#ifdef USEMPI
        MPITask.resize(NumberofDescendants);
        removalhaloindex.reserve(d.removalhaloindex.size());
        removalhalotemporalindex.reserve(d.removalhaloindex.size());
#endif
        for (int i=0;i<NumberofDescendants;i++) {
            haloindex[i]=d.haloindex[i];
            halotemporalindex[i]=d.halotemporalindex[i];
            Merit[i]=d.Merit[i];
            deltat[i]=d.deltat[i];
            descentype[i]=d.descentype[i];
#ifdef USEMPI
            MPITask[i]=d.MPITask[i];
#endif
        }
#ifdef USEMPI
        for (int i=0;i<d.removalhaloindex.size();i++) {
            removalhaloindex[i]=d.removalhaloindex[i];
            removalhalotemporalindex[i]=d.removalhalotemporalindex[i];
        }
#endif
    }
    ///Determine optimal descendent using a temporally weighted merit, and set it to position 0
    ///start with just maximum merit, ignoring when this was found
    ///otherwise use the reference time passed
    ///the generalized merit = Merit/(deltat)
    ///also, optimal descendents should be close to being a primary descendant, having a lower descentype value
    void OptimalTemporalMerit(int iopttemporalmerittype=GENERALIZEDMERITTIME, Int_t itimeref=0){
        int imax=0;
        Double_t generalizedmerit=Merit[0]/pow((Double_t)deltat[0],ALPHADELTAT), newgenmerit;
        int curdescentype=descentype[0];
        long unsigned optimalhaloindex;
        int unsigned optimalhalotemporalindex;

        for (int i=1;i<NumberofDescendants;i++) {
            newgenmerit=Merit[i]/pow((Double_t)deltat[i],ALPHADELTAT);
            //if just optimising generalized temporal merit
            if (iopttemporalmerittype==GENERALIZEDMERITTIME && newgenmerit>generalizedmerit) {
                generalizedmerit=newgenmerit;curdescentype=descentype[i];imax=i;
            }
            //if optimising for best merit and best merit ranking (that is how close the object is to being the primary progenitor)
            else if (iopttemporalmerittype==GENERALIZEDMERITTIMEPROGEN && ((descentype[i]<curdescentype && newgenmerit>=generalizedmerit*0.25)||(newgenmerit>generalizedmerit))) {
                generalizedmerit=newgenmerit;curdescentype=descentype[i];imax=i;
            }
        }
        optimalhaloindex=haloindex[imax];
        optimalhalotemporalindex=halotemporalindex[imax];
        //if use mpi then also possible that optimal halo found on multiple mpi tasks
        //but the lower task will be the one that needs to have the best halo so go over the loop and
        //find the halo with the lowest task number of the match
#ifdef USEMPI
        int imaxtask=MPITask[imax];
        for (int i=0;i<NumberofDescendants;i++) {
            if (optimalhalotemporalindex==halotemporalindex[i] && optimalhaloindex==haloindex[i] && MPITask[i]<imaxtask) {imax=i;}
        }
#endif
        if (imax>0) {
            long unsigned hid, htid;
            Double_t merit;
            int dt;
            int dtype;
#ifdef USEMPI
            int mpitask;
#endif
            merit=Merit[0];
            hid=haloindex[0];
            htid=halotemporalindex[0];
            dt=deltat[0];
            dtype=descentype[0];
#ifdef USEMPI
            mpitask=MPITask[0];
#endif
            Merit[0]=Merit[imax];
            haloindex[0]=haloindex[imax];
            halotemporalindex[0]=halotemporalindex[imax];
            deltat[0]=deltat[imax];
            descentype[0]=descentype[imax];
#ifdef USEMPI
            MPITask[0]=MPITask[imax];
#endif
            Merit[imax]=merit;
            haloindex[imax]=hid;
            halotemporalindex[imax]=htid;
            deltat[imax]=dt;
            descentype[imax]=dtype;
#ifdef USEMPI
            MPITask[imax]=mpitask;
#endif
        }
    }
#ifdef USEMPI
    void Merge(int thistask, int &numdescen, long unsigned *hid,int unsigned *htid, float *m, int *dt, int *dtype, int *task) {
        for (Int_t i=0;i<numdescen;i++) if (task[i]!=thistask)
        {
            //first check to see if halo already present
            int k=0;
            while (k<NumberofDescendants && !(haloindex[k]==hid[i] && halotemporalindex[k]==htid[i])) k++;
            //if entry not found do nothing
            if (k<NumberofDescendants) continue;
            haloindex.push_back(hid[i]);
            halotemporalindex.push_back(htid[i]);
            Merit.push_back(m[i]);
            deltat.push_back(dt[i]);
            descentype.push_back(dtype[i]);
            MPITask.push_back(task[i]);
            NumberofDescendants++;
        }
    }
    //remove entries from the list
    void Removal(int &nremove, long unsigned *hid,int unsigned *htid) {
        for (Int_t i=0;i<nremove;i++)
        {
            //find entry
            int k=0;
            while (k<NumberofDescendants && !(haloindex[k]==hid[i] && halotemporalindex[k]==htid[i])) k++;
            //if entry not found do nothing
            if (k==NumberofDescendants) continue;
            haloindex.erase(haloindex.begin()+k);
            halotemporalindex.erase(halotemporalindex.begin()+k);
            Merit.erase(Merit.begin()+k);
            deltat.erase(deltat.begin()+k);
            descentype.erase(descentype.begin()+k);
            MPITask.erase(MPITask.begin()+k);
            NumberofDescendants--;
        }
    }
#endif

};

/*!
    Structure used to keep track of a structure's candidate preogenitors based on the descendant construction
*/
struct ProgenitorDataDescenBased
{
    ///structure type and number of links (progenitors)
    //@{
    int stype;
    int NumberofProgenitors;
    //@}

    ///store list of progenitors in the form of halo index and temporal index
    vector<long unsigned> haloindex;
    vector<int unsigned> halotemporalindex;
    ///store the merit value
    vector<float> Merit;
    //store the integer time diff to progenitor
    vector<int> deltat;
    //store index in progenitor's descedant based list
    vector<int> progenindex;
    //store descendant to progenitor ranking (dtoptype)
    vector<int> dtoptype;
#ifdef USEMPI
    //store local number of mpi items
    int nlocal;
    //store which task this halo progenitor is located on
    vector<int> MPITask;
    //store items removed
    vector<long unsigned> removalhaloindex;
    vector<int unsigned> removalhalotemporalindex;
#endif
    ProgenitorDataDescenBased(int reservesize=4){
        NumberofProgenitors=0;
        //reserve some space so mimimize number of reallocations
        haloindex.reserve(reservesize);
        halotemporalindex.reserve(reservesize);
        Merit.reserve(reservesize);
        deltat.reserve(reservesize);
        progenindex.reserve(reservesize);
        dtoptype.reserve(reservesize);
#ifdef USEMPI
        MPITask.reserve(reservesize);
        removalhaloindex.reserve(reservesize);
        removalhalotemporalindex.reserve(reservesize);
#endif
    }
    ~ProgenitorDataDescenBased(){
    }
    ProgenitorDataDescenBased &operator=(const ProgenitorDataDescenBased &d){
        NumberofProgenitors=d.NumberofProgenitors;
        haloindex.resize(NumberofProgenitors);
        halotemporalindex.resize(NumberofProgenitors);
        Merit.resize(NumberofProgenitors);
        progenindex.reserve(NumberofProgenitors);
        deltat.resize(NumberofProgenitors);
        dtoptype.reserve(NumberofProgenitors);
#ifdef USEMPI
        MPITask.resize(NumberofProgenitors);
        removalhaloindex.reserve(d.removalhaloindex.size());
        removalhalotemporalindex.reserve(d.removalhaloindex.size());
#endif
        for (int i=0;i<NumberofProgenitors;i++) {
            haloindex[i]=d.haloindex[i];
            halotemporalindex[i]=d.halotemporalindex[i];
            Merit[i]=d.Merit[i];
            deltat[i]=d.deltat[i];
            progenindex[i]=d.progenindex[i];
            dtoptype[i]=d.dtoptype[i];
#ifdef USEMPI
            MPITask[i]=d.MPITask[i];
#endif
        }
#ifdef USEMPI
        for (int i=0;i<d.removalhaloindex.size();i++) {
            removalhaloindex[i]=d.removalhaloindex[i];
            removalhalotemporalindex[i]=d.removalhalotemporalindex[i];
        }
#endif
    }
    ///Determine optimal descendent using a temporally weighted merit, and set it to position 0
    ///start with just maximum merit, ignoring when this was found
    ///otherwise use the reference time passed
    ///the generalized merit = Merit/(deltat)
    ///also, optimal descendents should be close to being a primary progenitor, having a lower descentype value
    void OptimalTemporalMerit(int iopttemporalmerittype=GENERALIZEDMERITTIME, Int_t itimeref=0){
        int imax=0;
        Double_t generalizedmerit=Merit[0]/pow((Double_t)deltat[0],ALPHADELTAT), newgenmerit;
        int curprogenindex=progenindex[0];
        long unsigned optimalhaloindex;
        int unsigned optimalhalotemporalindex;
        for (int i=1;i<NumberofProgenitors;i++) {
            newgenmerit=Merit[i]/pow((Double_t)deltat[i],ALPHADELTAT);
            //if just optimising generalized temporal merit
            if (iopttemporalmerittype==GENERALIZEDMERITTIME && newgenmerit>generalizedmerit) {
                generalizedmerit=newgenmerit;curprogenindex=progenindex[i];imax=i;
            }
            //if optimising for best merit and best merit ranking (that is how close the object is to being the primary progenitor)
            else if (iopttemporalmerittype==GENERALIZEDMERITTIMEPROGEN && ((progenindex[i]<curprogenindex && newgenmerit>=generalizedmerit*0.25)||(newgenmerit>generalizedmerit))) {
                generalizedmerit=newgenmerit;curprogenindex=progenindex[i];imax=i;
            }
        }
        optimalhaloindex=haloindex[imax];
        optimalhalotemporalindex=halotemporalindex[imax];
        //if use mpi then also possible that optimal halo found on multiple mpi tasks
        //but the lower task will be the one that needs to have the best halo so go over the loop and
        //find the halo with the lowest task number of the best generalized merit
#ifdef USEMPI
        int imaxtask=MPITask[imax];
        for (int i=0;i<NumberofProgenitors;i++) {
            if (optimalhalotemporalindex==halotemporalindex[i] && optimalhaloindex==haloindex[i] && MPITask[i]<imaxtask) {imax=i;}
        }
#endif
        if (imax>0) {
            long unsigned hid, htid;
            Double_t merit;
            int dt;
            int pindex;
            int dtoptypeval;
#ifdef USEMPI
            int mpitask;
#endif
            merit=Merit[0];
            hid=haloindex[0];
            htid=halotemporalindex[0];
            dt=deltat[0];
            pindex=progenindex[0];
            dtoptypeval=dtoptype[0];
#ifdef USEMPI
            mpitask=MPITask[0];
#endif
            Merit[0]=Merit[imax];
            haloindex[0]=haloindex[imax];
            halotemporalindex[0]=halotemporalindex[imax];
            deltat[0]=deltat[imax];
            progenindex[0]=progenindex[imax];
#ifdef USEMPI
            MPITask[0]=MPITask[imax];
#endif
            Merit[imax]=merit;
            haloindex[imax]=hid;
            halotemporalindex[imax]=htid;
            deltat[imax]=dt;
            progenindex[imax]=pindex;
            dtoptype[imax]=dtoptypeval;
#ifdef USEMPI
            MPITask[imax]=mpitask;
#endif
        }
    }
#ifdef USEMPI
    void Merge(int thistask, int &numprogen, long unsigned *hid,int unsigned *htid, float *m, int *dt, int *pindex, int* dtop, int *task) {
        for (Int_t i=0;i<numprogen;i++) if (task[i]!=thistask)
        {
            //first check to see if halo already present
            int k=0;
            while (k<NumberofProgenitors && !(haloindex[k]==hid[i] && halotemporalindex[k]==htid[i])) k++;
            //if entry not found do nothing
            if (k<NumberofProgenitors) continue;
            haloindex.push_back(hid[i]);
            halotemporalindex.push_back(htid[i]);
            Merit.push_back(m[i]);
            deltat.push_back(dt[i]);
            progenindex.push_back(pindex[i]);
            dtoptype.push_back(dtop[i]);
            MPITask.push_back(task[i]);
            NumberofProgenitors++;
        }
    }
    //remove entries from the list
    void Removal(int &nremove, long unsigned *hid,int unsigned *htid) {
        for (Int_t i=0;i<nremove;i++)
        {
            //find entry
            int k=0;
            while (k<NumberofProgenitors && !(haloindex[k]==hid[i] && halotemporalindex[k]==htid[i])) k++;
            //if entry not found do nothing
            if (k==NumberofProgenitors) continue;
            haloindex.erase(haloindex.begin()+k);
            halotemporalindex.erase(halotemporalindex.begin()+k);
            Merit.erase(Merit.begin()+k);
            deltat.erase(deltat.begin()+k);
            progenindex.erase(progenindex.begin()+k);
            dtoptype.erase(dtoptype.begin()+k);
            MPITask.erase(MPITask.begin()+k);
            NumberofProgenitors--;
        }
    }
#endif

};


///struct to store the file type names produced by velociraptor
struct VELOCIraptorFileTypeNames {
    vector<string> ftypename;
    VELOCIraptorFileTypeNames(){
        ftypename.push_back("catalog_groups");
        ftypename.push_back("catalog_particles");
        ftypename.push_back("catalog_particles.unbound");
        ftypename.push_back("sublevels.catalog_groups");
        ftypename.push_back("sublevels.catalog_particles");
        ftypename.push_back("sublevels.catalog_particles.unbound");
        ftypename.push_back("catalog_parttypes");
        ftypename.push_back("catalog_parttypes.unbound");
        ftypename.push_back("sublevels.catalog_parttypes");
        ftypename.push_back("sublevels.catalog_parttypes.unbound");
        ftypename.push_back("properties");
        ftypename.push_back("sublevels.properties");
    }
    void UpdateName(string &infile, string farray[], int k,int impi=0){
        std::ostringstream oss;
        oss<<"."<<k;
        for (int i=0;i<ftypename.size();i++) {farray[i]=infile;farray[i]+=".";farray[i]+=ftypename[i];}
        if (impi) for (int i=0;i<ftypename.size();i++) {farray[i]+=oss.str();}
    }
};

#ifdef USEHDF
///hdf structure to store the names of datasets in catalog output
struct HDFCatalogNames {
    ///store names of catalog group files
    vector<H5std_string> group;
    //store the data type
    vector<PredType> groupdatatype;

    ///store the names of catalog particle files
    vector<H5std_string> part;
    //store the data type
    vector<PredType> partdatatype;

    ///store the names of catalog particle files
    vector<H5std_string> types;
    //store the data type
    vector<PredType> typesdatatype;

    ///store the names of catalog particle files
    vector<H5std_string> hierarchy;
    //store the data type
    vector<PredType> hierarchydatatype;

    ///store the names of catalog particle files
    vector<H5std_string> properties;
    //store the data type
    vector<PredType> propertiesdatatype;

    HDFCatalogNames(){
        group.push_back("File_id");
        group.push_back("Num_of_files");
        group.push_back("Num_of_groups");
        group.push_back("Total_num_of_groups");
        group.push_back("Group_Size");
        group.push_back("Offset");
        group.push_back("Offset_unbound");
        groupdatatype.push_back(PredType::STD_I32LE);
        groupdatatype.push_back(PredType::STD_I32LE);
        groupdatatype.push_back(PredType::STD_U64LE);
        groupdatatype.push_back(PredType::STD_U64LE);
        groupdatatype.push_back(PredType::STD_U32LE);
        groupdatatype.push_back(PredType::STD_U64LE);
        groupdatatype.push_back(PredType::STD_U64LE);

        part.push_back("File_id");
        part.push_back("Num_of_files");
        part.push_back("Num_of_particles_in_groups");
        part.push_back("Total_num_of_particles_in_all_groups");
        part.push_back("Particle_IDs");
        partdatatype.push_back(PredType::STD_I32LE);
        partdatatype.push_back(PredType::STD_I32LE);
        partdatatype.push_back(PredType::STD_U64LE);
        partdatatype.push_back(PredType::STD_U64LE);
        partdatatype.push_back(PredType::STD_I64LE);

        types.push_back("File_id");
        types.push_back("Num_of_files");
        types.push_back("Num_of_particles_in_groups");
        types.push_back("Total_num_of_particles_in_all_groups");
        types.push_back("Particle_types");
        typesdatatype.push_back(PredType::STD_I32LE);
        typesdatatype.push_back(PredType::STD_I32LE);
        typesdatatype.push_back(PredType::STD_U64LE);
        typesdatatype.push_back(PredType::STD_U64LE);
        typesdatatype.push_back(PredType::STD_U16LE);


        hierarchy.push_back("File_id");
        hierarchy.push_back("Num_of_files");
        hierarchy.push_back("Num_of_groups");
        hierarchy.push_back("Total_num_of_groups");
        hierarchy.push_back("Number_of_substructures_in_halo");
        hierarchy.push_back("Parent_halo_ID");
        hierarchydatatype.push_back(PredType::STD_I32LE);
        hierarchydatatype.push_back(PredType::STD_I32LE);
        hierarchydatatype.push_back(PredType::STD_U64LE);
        hierarchydatatype.push_back(PredType::STD_U64LE);
        hierarchydatatype.push_back(PredType::STD_U32LE);
        hierarchydatatype.push_back(PredType::STD_I64LE);

        properties.push_back("File_id");
        properties.push_back("Num_of_files");
        properties.push_back("Num_of_groups");
        properties.push_back("Total_num_of_groups");
        properties.push_back("Xc");
        properties.push_back("Yc");
        properties.push_back("Zc");
        properties.push_back("VXc");
        properties.push_back("VYc");
        properties.push_back("VZc");
        properties.push_back("Rmax");
        properties.push_back("Vmax");
        propertiesdatatype.push_back(PredType::STD_I32LE);
        propertiesdatatype.push_back(PredType::STD_I32LE);
        propertiesdatatype.push_back(PredType::STD_U64LE);
        propertiesdatatype.push_back(PredType::STD_U64LE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
        propertiesdatatype.push_back(PredType::NATIVE_DOUBLE);
    }
};
#endif

#endif
