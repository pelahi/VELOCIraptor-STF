/*! \file buildandsortarrays.cxx
 *  \brief this file contains routines that build arrays used to sort/access the particle data local to the MPI domain
 */

#include "stf.h"

/// \name Simple group id based array building and group id reordering routines
//@{

///build group size array
Int_t *BuildNumInGroup(const Int_t nbodies, const Int_t numgroups, Int_t *pfof){
    Int_t *numingroup=new Int_t[numgroups+1];
    for (Int_t i=0;i<=numgroups;i++) numingroup[i]=0;
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) numingroup[pfof[i]]++;
    return numingroup;
}
///build group size array for specific type
Int_t *BuildNumInGroupTyped(const Int_t nbodies, const Int_t numgroups, Int_t *pfof, Particle *P, int type){
    Int_t *numingroup=new Int_t[numgroups+1];
    for (Int_t i=0;i<=numgroups;i++) numingroup[i]=0;
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0 && P[i].GetType()==type) numingroup[pfof[i]]++;
    return numingroup;
}

///build the group particle index list (assumes particles are in ID order)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof){
    Int_t **pglist=new Int_t*[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++) {pglist[i]=new Int_t[numingroup[i]+1*(numingroup[i]==0)];numingroup[i]=0;}
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) pglist[pfof[i]][numingroup[pfof[i]]++]=i;
    return pglist;
}
///build the group particle index list for particles of a specific type (assumes particles are in ID order)
Int_t **BuildPGListTyped(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *P, int type){
    Int_t **pglist=new Int_t*[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++) {pglist[i]=new Int_t[numingroup[i]+1*(numingroup[i]==0)];numingroup[i]=0;}
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0 && P[i].GetType()==type) pglist[pfof[i]][numingroup[pfof[i]]++]=i;
    return pglist;
}
///build the group particle index list (doesn't assume particles are in ID order and stores index of particle)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Particle *Part){
    Int_t **pglist=new Int_t*[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++) {pglist[i]=new Int_t[numingroup[i]];numingroup[i]=0;}
    for (Int_t i=0;i<nbodies;i++) {
        int pid=Part[i].GetID();
        if (pfof[pid]>0) pglist[pfof[pid]][numingroup[pfof[pid]]++]=i;
    }
    return pglist;
}
///build the group particle index list (doesn't assumes particles are in ID order)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof, Int_t *ids){
    Int_t **pglist=new Int_t*[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++) {pglist[i]=new Int_t[numingroup[i]];numingroup[i]=0;}
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) pglist[pfof[i]][numingroup[pfof[i]]++]=ids[i];
    return pglist;
}
///build the Head array which points to the head of the group a particle belongs to
Int_tree_t *BuildHeadArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *Head=new Int_tree_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) Head[i]=i;
    for (Int_t i=1;i<=numgroups;i++)
        for (Int_t j=1;j<numingroup[i];j++) Head[pglist[i][j]]=Head[pglist[i][0]];
    return Head;
}
///build the Next array which points to the next particle in the group
Int_tree_t *BuildNextArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *Next=new Int_tree_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) Next[i]=-1;
    for (Int_t i=1;i<=numgroups;i++)
        for (Int_t j=0;j<numingroup[i]-1;j++) Next[pglist[i][j]]=pglist[i][j+1];
    return Next;
}
///build the Len array which stores the length of the group a particle belongs to
Int_tree_t *BuildLenArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *Len=new Int_tree_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) Len[i]=0;
    for (Int_t i=1;i<=numgroups;i++)
        for (Int_t j=0;j<numingroup[i];j++) Len[pglist[i][j]]=numingroup[i];
    return Len;
}
///build the GroupTail array which stores the Tail of a group
Int_tree_t *BuildGroupTailArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_tree_t *GTail=new Int_tree_t[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++)
        GTail[i]=pglist[i][numingroup[i]-1];
    return GTail;
}
///build the group particle arrays need for unbinding procedure
Particle **BuildPartList(Int_t numgroups, Int_t *numingroup, Int_t **pglist, Particle* Part){
    Particle **gPart=new Particle*[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++) {
        gPart[i]=new Particle[numingroup[i]];
        for (Int_t j=0;j<numingroup[i];j++) gPart[i][j]=Part[pglist[i][j]];
    }
    return gPart;
}
///build a particle list subset using array of indices
Particle *BuildPart(Int_t numingroup, Int_t *pglist, Particle* Part){
    Particle *gPart=new Particle[numingroup+1];
    for (Int_t j=0;j<numingroup;j++) gPart[j]=Part[pglist[j]];
    return gPart;
}

///sort particles according to some quantity which is stored in particle type and build an array for a sorted particle list
///remember this reorders the particle array!
Int_t *BuildNoffset(const Int_t nbodies, Particle *Part, Int_t numgroups,Int_t *numingroup, Int_t *sortval, Int_t ioffset) {
    Int_t *noffset=new Int_t[numgroups+1];
    Int_t *storetype=new Int_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) {
        storetype[Part[i].GetID()]=Part[i].GetType();
        if (sortval[Part[i].GetID()]>ioffset) Part[i].SetType(sortval[Part[i].GetID()]);
        else Part[i].SetType(nbodies+1);//here move all particles not in groups to the back of the particle array
    }
    qsort(Part, nbodies, sizeof(Particle), TypeCompare);
    if (numgroups > 1) noffset[0]=noffset[1]=0;
    for (Int_t i=2;i<=numgroups;i++) noffset[i]=noffset[i-1]+numingroup[i-1];
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
    Int_t groupid,size;
    for (Int_t i = 1; i <=numgroups; i++) if (numingroup[i]>0) pq->Push(i, numingroup[i]);
    for (Int_t i = 1; i<=newnumgroups; i++) {
        groupid=pq->TopQueue();size=pq->TopPriority();pq->Pop();
        for (Int_t j=0;j<size;j++) pfof[pglist[groupid][j]]=i;
    }
    delete pq;
}
void ReorderGroupIDs(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Particle *Partsubset)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Int_t groupid,size;
    for (Int_t i = 1; i <=numgroups; i++) if (numingroup[i]>0) pq->Push(i, numingroup[i]);
    for (Int_t i = 1; i<=newnumgroups; i++) {
        groupid=pq->TopQueue();size=pq->TopPriority();pq->Pop();
        for (Int_t j=0;j<size;j++) pfof[Partsubset[pglist[groupid][j]].GetID()]=i;
    }
    delete pq;
}

///similar to \ref ReorderGroupIDs but weight by value
void ReorderGroupIDsbyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Int_t groupid;
    for (Int_t i = 1; i <=numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    for (Int_t i = 1; i<=newnumgroups; i++) {
        groupid=pq->TopQueue();pq->Pop();
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
    }
    delete pq;
}
///similar to \ref ReorderGroupIDsbyValue but also reorder associated group data
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Int_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Int_t *gtemp=new Int_t[numgroups+1];
    Int_t groupid;
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    for (Int_t i = 1; i <= newnumgroups; i++) {
        groupid=pq->TopQueue();pq->Pop();
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
        gdata[i]=gtemp[groupid];
    }
    delete pq;
    delete[] gtemp;
}
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Double_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Double_t *gtemp=new Double_t[numgroups+1];
    Int_t groupid;
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    for (Int_t i = 1; i <= newnumgroups; i++) {
        groupid=pq->TopQueue();pq->Pop();
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
        gdata[i]=gtemp[groupid];
    }
    delete pq;
    delete[] gtemp;
}
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Double_t *value, Int_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Double_t *gtemp=new Double_t[numgroups+1];
    Int_t groupid;
    Int_t count=0;
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    for (Int_t i = 1; i <= newnumgroups; i++) {
        groupid=pq->TopQueue();pq->Pop();
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
        gdata[i]=gtemp[groupid];
    }
    delete pq;
    delete[] gtemp;
}
void ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Double_t *value, Double_t *gdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    Double_t *gtemp=new Double_t[numgroups+1];
    Int_t groupid;
    Int_t count=0;
    for (Int_t i = 1; i <= numgroups; i++) gtemp[i]=gdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    for (Int_t i = 1; i <= newnumgroups; i++) {
        groupid=pq->TopQueue();pq->Pop();
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
        gdata[i]=gtemp[groupid];
    }
    delete pq;
    delete[] gtemp;
}
///similar to \ref ReorderGroupIDsbyValue but also reorder associated property data
void ReorderGroupIDsAndHaloDatabyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, PropData *pdata)
{
    PriorityQueue *pq=new PriorityQueue(newnumgroups);
    PropData *ptemp=new PropData[numgroups+1];
    Int_t groupid;
    for (Int_t i = 1; i <= numgroups; i++) ptemp[i]=pdata[i];
    for (Int_t i = 1; i <= numgroups; i++) if (numingroup[i]>0) pq->Push(i, value[i]);
    for (Int_t i = 1; i <= newnumgroups; i++) {
        groupid=pq->TopQueue();pq->Pop();
        for (Int_t j=0;j<numingroup[groupid];j++) pfof[pglist[groupid][j]]=i;
        pdata[i]=ptemp[groupid];
    }
    delete pq;
    delete[] ptemp;
}
//@}
