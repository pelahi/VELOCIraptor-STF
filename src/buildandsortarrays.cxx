/*! \file buildandsortarrays.cxx
 *  \brief this file contains routines that build arrays used to sort/access the particle data local to the MPI domain
 */

#include "stf.h"

#include <parallel/algorithm>

/// \name Simple group id based array building and group id reordering routines
//@{

///build group size array
Int_t *BuildNumInGroup(const Int_t nbodies, const Int_t numgroups, Int_t *pfof){
    Int_t *numingroup=new Int_t[numgroups+1];
    for (Int_t i=0;i<=numgroups;i++) numingroup[i]=0;
#ifdef USEOPENMP
    //guarded by if(!omp_in_parallel()): if this is already called from within an
    //active parallel region, this whole block runs as a team of one (no new threads,
    //no oversubscription) and just falls through to the plain sequential accumulation.
    #pragma omp parallel if(!omp_in_parallel())
    {
        Int_t *localcount = new Int_t[numgroups+1]();
        #pragma omp for schedule(static) nowait
        for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) localcount[pfof[i]]++;
        //per-thread histograms merged here instead of atomics on numingroup[pfof[i]]
        //directly, since a handful of huge groups would otherwise force heavy
        //cross-thread contention on the same few buckets.
        #pragma omp critical
        {
            for (Int_t i=0;i<=numgroups;i++) numingroup[i]+=localcount[i];
        }
        delete[] localcount;
    }
#else
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) numingroup[pfof[i]]++;
#endif
    return numingroup;
}
///build group size array for specific type
Int_t *BuildNumInGroupTyped(const Int_t nbodies, const Int_t numgroups, Int_t *pfof, Particle *P, int type){
    Int_t *numingroup=new Int_t[numgroups+1];
    for (Int_t i=0;i<=numgroups;i++) numingroup[i]=0;
#ifdef USEOPENMP
    #pragma omp parallel if(!omp_in_parallel())
    {
        Int_t *localcount = new Int_t[numgroups+1]();
        #pragma omp for schedule(static) nowait
        for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0 && P[i].GetType()==type) localcount[pfof[i]]++;
        #pragma omp critical
        {
            for (Int_t i=0;i<=numgroups;i++) numingroup[i]+=localcount[i];
        }
        delete[] localcount;
    }
#else
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0 && P[i].GetType()==type) numingroup[pfof[i]]++;
#endif
    return numingroup;
}

///Order-preserving parallel scatter into a pre-sized pglist: for i in [0,nbodies), if
///getpid(i) names a valid group, writes getval(i) into pglist[pid] at the position it
///would have landed in under the original strictly-sequential loop (i=0,1,2,...), i.e.
///with the same left-to-right particle order within each group that the single-threaded
///version produced. This matters because an atomic-capture "first thread there wins the
///next slot" scatter scrambles that per-group order based on thread scheduling, which
///was observed to perturb an order-sensitive downstream fit (NFW concentration).
///Achieved via a counting pass (schedule(static), so thread t always owns a fixed,
///contiguous, increasing-index chunk) followed by a per-group prefix sum across threads
///to hand each thread a deterministic non-overlapping write range, then a second pass
///that writes using those ranges -- no atomics needed since ranges never overlap.
template<typename GetPidFn, typename GetValFn>
static void ScatterIntoPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup,
                               Int_t **pglist, GetPidFn getpid, GetValFn getval)
{
#ifdef USEOPENMP
    int nthreads=1;
    if (!omp_in_parallel()) {
        #pragma omp parallel
        {
            #pragma omp single
            nthreads = omp_get_num_threads();
        }
    }
    vector<vector<Int_t>> localcount(nthreads, vector<Int_t>(numgroups+1,0));
    #pragma omp parallel if(!omp_in_parallel()) num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        #pragma omp for schedule(static)
        for (Int_t i=0;i<nbodies;i++) {
            Int_t pid = getpid(i);
            if (pid==0 || numingroup[pid]<0) continue;
            localcount[tid][pid]++;
        }
    }
    vector<vector<Int_t>> offset(nthreads, vector<Int_t>(numgroups+1,0));
    for (Int_t pid=1; pid<=numgroups; pid++) {
        Int_t running=0;
        for (int t=0;t<nthreads;t++) { offset[t][pid]=running; running+=localcount[t][pid]; }
    }
    #pragma omp parallel if(!omp_in_parallel()) num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        vector<Int_t> cursor = offset[tid];
        #pragma omp for schedule(static)
        for (Int_t i=0;i<nbodies;i++) {
            Int_t pid = getpid(i);
            if (pid==0 || numingroup[pid]<0) continue;
            pglist[pid][cursor[pid]++]=getval(i);
        }
    }
#else
    vector<Int_t> cursor(numgroups+1,0);
    for (Int_t i=0;i<nbodies;i++) {
        Int_t pid = getpid(i);
        if (pid==0 || numingroup[pid]<0) continue;
        pglist[pid][cursor[pid]++]=getval(i);
    }
#endif
}

///build the group particle index list (assumes particles are in ID order)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof){
    Int_t **pglist=new Int_t*[numgroups+1];
    pglist[0]=NULL;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        pglist[i] = NULL;
        if (numingroup[i]<=0) continue;
        pglist[i]=new Int_t[numingroup[i]];
    }
    ScatterIntoPGList(nbodies, numgroups, numingroup, pglist,
        [pfof](Int_t i){ return pfof[i]; },
        [](Int_t i){ return i; });
    return pglist;
}
///build the group particle index list for particles of a specific type (assumes particles are in ID order)
Int_t **BuildPGListTyped(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *P, int type){
    Int_t **pglist=new Int_t*[numgroups+1];
    pglist[0]=NULL;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        pglist[i] = NULL;
        if (numingroup[i]<=0) continue;
        pglist[i]=new Int_t[numingroup[i]];
    }
    ScatterIntoPGList(nbodies, numgroups, numingroup, pglist,
        [pfof,P,type](Int_t i){ return (P[i].GetType()==type) ? pfof[i] : (Int_t)0; },
        [](Int_t i){ return i; });
    return pglist;
}
///build the group particle index list (doesn't assume particles are in ID order and stores index of particle)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *Part){
    Int_t **pglist=new Int_t*[numgroups+1];
    pglist[0]=NULL;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        pglist[i] = NULL;
        if (numingroup[i]<=0) continue;
        pglist[i]=new Int_t[numingroup[i]];
    }
    ScatterIntoPGList(nbodies, numgroups, numingroup, pglist,
        [pfof,Part](Int_t i){ return pfof[Part[i].GetID()]; },
        [](Int_t i){ return i; });
    return pglist;
}
///build the group particle index list (doesn't assumes particles are in ID order)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Int_t *ids){
    Int_t **pglist=new Int_t*[numgroups+1];
    pglist[0]=NULL;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        pglist[i] = NULL;
        if (numingroup[i]<=0) continue;
        pglist[i]=new Int_t[numingroup[i]];
    }
    ScatterIntoPGList(nbodies, numgroups, numingroup, pglist,
        [pfof](Int_t i){ return pfof[i]; },
        [ids](Int_t i){ return ids[i]; });
    return pglist;
}
///build the Head array which points to the head of the group a particle belongs to
Int_tree_t *BuildHeadArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *Head=new Int_tree_t[nbodies];
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=0;i<nbodies;i++) Head[i]=i;
    //each group only ever writes to indices named in its own pglist[i], disjoint from
    //every other group's, so this is safe to parallelize directly across groups.
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        for (Int_t j=1;j<numingroup[i];j++) Head[pglist[i][j]]=Head[pglist[i][0]];
    }
    return Head;
}
///build the Next array which points to the next particle in the group
Int_tree_t *BuildNextArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *Next=new Int_tree_t[nbodies];
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=0;i<nbodies;i++) Next[i]=-1;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        for (Int_t j=0;j<numingroup[i]-1;j++) Next[pglist[i][j]]=pglist[i][j+1];
    }
    return Next;
}
///build the Len array which stores the length of the group a particle belongs to
Int_tree_t *BuildLenArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *Len=new Int_tree_t[nbodies];
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=0;i<nbodies;i++) Len[i]=0;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        for (Int_t j=0;j<numingroup[i];j++) Len[pglist[i][j]]=numingroup[i];
    }
    return Len;
}
///build the GroupTail array which stores the Tail of a group
Int_tree_t *BuildGroupTailArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *GTail=new Int_tree_t[numgroups+1];
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=1;i<=numgroups;i++) {
        GTail[i]=pglist[i][numingroup[i]-1];
    }
    return GTail;
}
///build the group particle arrays need for unbinding procedure
Particle **BuildPartList(Int_t numgroups, Int_t *numingroup, Int_t **pglist, Particle* Part,
    bool ikeepextrainfo)
{
    Particle **gPart=new Particle*[numgroups+1];
    gPart[0] = NULL;
    //each group's gPart[i] is independent of every other group's, so this is safe to
    //parallelize directly across groups with no atomics needed; guarded the same way
    //as BuildNumInGroup/BuildPGList against an already-active outer parallel region.
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (auto i=1;i<=numgroups;i++) {
        gPart[i] = NULL;
        if (numingroup[i]<=0) continue;
        gPart[i]=new Particle[numingroup[i]];
        for (auto j=0;j<numingroup[i];j++) {
            gPart[i][j]=Part[pglist[i][j]];
#if defined(GASON) || defined(STARON) || defined(BHON) || defined(EXTRADMON)
            if (ikeepextrainfo) continue;
#endif
#ifdef GASON
            if (gPart[i][j].HasHydroProperties()) gPart[i][j].SetHydroProperties();
#endif
#ifdef STARON
            if (gPart[i][j].HasStarProperties()) gPart[i][j].SetStarProperties();
#endif
#ifdef BHON
            if (gPart[i][j].HasBHProperties()) gPart[i][j].SetBHProperties();
#endif
#ifdef EXTRADMON
            if (gPart[i][j].HasExtraDMProperties()) gPart[i][j].SetExtraDMProperties();
#endif
        }
    }
    return gPart;
}
///build a particle list subset using array of indices
Particle *BuildPart(Int_t numingroup, Int_t *pglist, Particle* Part,
    bool ikeepextrainfo)
{
    Particle *gPart=new Particle[numingroup+1];
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (auto j=0;j<numingroup;j++) {
        gPart[j]=Part[pglist[j]];
#if defined(GASON) || defined(STARON) || defined(BHON) || defined(EXTRADMON)
            if (ikeepextrainfo) continue;
#endif
#ifdef GASON
        if (gPart[j].HasHydroProperties()) gPart[j].SetHydroProperties();
#endif
#ifdef STARON
        if (gPart[j].HasStarProperties()) gPart[j].SetStarProperties();
#endif
#ifdef BHON
        if (gPart[j].HasBHProperties()) gPart[j].SetBHProperties();
#endif
#ifdef EXTRADMON
        if (gPart[j].HasExtraDMProperties()) gPart[j].SetExtraDMProperties();
#endif
    }
    return gPart;
}

///sort particles according to some quantity which is stored in particle type and build an array for a sorted particle list
///remember this reorders the particle array!
Int_t *BuildNoffset(const Int_t nbodies, Particle *Part, Int_t numgroups,Int_t *numingroup, Int_t *sortval, Int_t ioffset) {
    Int_t *noffset=new Int_t[numgroups+1];
    Int_t *storetype=new Int_t[nbodies];
    //GetID() ranges over [0,nbodies) as a permutation (relied on below to index
    //storetype/Part back and forth), so each iteration here touches disjoint elements
    //of both arrays and is safe to run in parallel.
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=0;i<nbodies;i++) {
        storetype[Part[i].GetID()]=Part[i].GetType();
        if (sortval[Part[i].GetID()]>ioffset) Part[i].SetType(sortval[Part[i].GetID()]);
        else Part[i].SetType(nbodies+1);//here move all particles not in groups to the back of the particle array
    }
    //was a single-threaded qsort; __gnu_parallel::stable_sort with an inlined ascending
    //comparator (matching TypeCompare's ordering) is dramatically faster for large
    //nbodies, same fix as applied to the id-order sort in search.cxx. Must be the
    //*stable* variant: Type here is set to a group id shared by every particle in
    //that group (a handful of distinct values, massively duplicated), so an unstable
    //sort would scramble each group's internal particle order.
#ifdef USEOPENMP
    __gnu_parallel::stable_sort(Part, Part+nbodies, [](const Particle &a, const Particle &b){ return a.GetType() < b.GetType(); });
#else
    qsort(Part, nbodies, sizeof(Particle), TypeCompare);
#endif
    if (numgroups >= 1) noffset[0]=noffset[1]=0;
    for (Int_t i=2;i<=numgroups;i++) noffset[i]=noffset[i-1]+numingroup[i-1];
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(static)
#endif
    for (Int_t i=0;i<nbodies;i++) Part[i].SetType(storetype[Part[i].GetID()]);
    delete[] storetype;
    return noffset;
}

///reorder groups from largest to smallest
///\todo must alter so that after pfof is reorderd, so is numingroup array and pglist so that do not have to reconstruct this list
///after reordering if numgroups==newnumgroups (ie, list has not shrunk)
void ReorderGroupIDs(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    for (Int_t i = 1; i <=numgroups; i++) if (numingroup[i]>0) pq->Push(i, numingroup[i]);
    //popping the PQ must stay sequential (priority order determines the new group id),
    //but it's only O(numgroups); capture (groupid,size) per new id here first, then do
    //the actual nbodies-scale pfof writes below in parallel -- each old group's
    //particles are disjoint from every other group's, so that part is safe to split
    //across threads regardless of pop order.
    Int_t *grouporder = new Int_t[newnumgroups+1];
    Int_t *groupsize = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i<=newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();groupsize[i]=pq->TopPriority();pq->Pop();
    }
    delete pq;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i<=newnumgroups; i++) {
        Int_t groupid=grouporder[i], size=groupsize[i];
        for (Int_t j=0;j<size;j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
    delete[] groupsize;
}
void ReorderGroupIDs(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Particle *Partsubset)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    for (Int_t i = 1; i <=numgroups; i++) if (numingroup[i]>0) pq->Push(i, numingroup[i]);
    Int_t *grouporder = new Int_t[newnumgroups+1];
    Int_t *groupsize = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i<=newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();groupsize[i]=pq->TopPriority();pq->Pop();
    }
    delete pq;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i<=newnumgroups; i++) {
        Int_t groupid=grouporder[i], size=groupsize[i];
        for (Int_t j=0;j<size;j++) pfof[Partsubset[pglist[groupid][j]].GetID()]=i;
    }
    delete[] grouporder;
    delete[] groupsize;
}

///similar to \ref ReorderGroupIDs but weight by value
void ReorderGroupIDsbyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    for (Int_t i = 1; i <=numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    Int_t *grouporder = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i<=newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();pq->Pop();
    }
    delete pq;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i<=newnumgroups; i++) {
        Int_t groupid=grouporder[i];
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
}
///similar to \ref ReorderGroupIDsbyValue but also reorder associated group data
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Int_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Int_t *gtemp=new Int_t[numgroups+1];
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    Int_t *grouporder = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i <= newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();pq->Pop();
        gdata[i]=gtemp[grouporder[i]];
    }
    delete pq;
    delete[] gtemp;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i <= newnumgroups; i++) {
        Int_t groupid=grouporder[i];
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
}
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Double_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Double_t *gtemp=new Double_t[numgroups+1];
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    Int_t *grouporder = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i <= newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();pq->Pop();
        gdata[i]=gtemp[grouporder[i]];
    }
    delete pq;
    delete[] gtemp;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i <= newnumgroups; i++) {
        Int_t groupid=grouporder[i];
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
}
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Double_t *value, Int_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Double_t *gtemp=new Double_t[numgroups+1];
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    //see ReorderGroupIDs: PQ popping stays sequential (cheap, O(numgroups)); the
    //nbodies-scale pfof writes are pulled out below to run in parallel.
    Int_t *grouporder = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i <= newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();pq->Pop();
        gdata[i]=gtemp[grouporder[i]];
    }
    delete pq;
    delete[] gtemp;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i <= newnumgroups; i++) {
        Int_t groupid=grouporder[i];
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
}
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Double_t *value, Double_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Double_t *gtemp=new Double_t[numgroups+1];
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    Int_t *grouporder = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i <= newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();pq->Pop();
        gdata[i]=gtemp[grouporder[i]];
    }
    delete pq;
    delete[] gtemp;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i <= newnumgroups; i++) {
        Int_t groupid=grouporder[i];
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
}
///similar to \ref ReorderGroupIDsbyValue but also reorder associated property data
void ReorderGroupIDsAndHaloDatabyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, PropData *pdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    PropData *ptemp=new PropData[numgroups+1];
    for (Int_t i = 1; i <= numgroups; i++) ptemp[i]=pdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    Int_t *grouporder = new Int_t[newnumgroups+1];
    for (Int_t i = 1; i <= newnumgroups; i++) {
        grouporder[i]=pq->TopQueue();pq->Pop();
        pdata[i]=ptemp[grouporder[i]];
    }
    delete pq;
    delete[] ptemp;
#ifdef USEOPENMP
    #pragma omp parallel for if(!omp_in_parallel()) schedule(dynamic)
#endif
    for (Int_t i = 1; i <= newnumgroups; i++) {
        Int_t groupid=grouporder[i];
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete[] grouporder;
}
//@}
