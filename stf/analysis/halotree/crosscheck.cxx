
/*! \file crosscheck.cxx
 *  \brief this file contains routines for determining and calculating cross matches
 */

#include "halomergertree.h"


/// \name cross matching routines
//@{

/// Determine initial progenitor list based on just merit
/// does not guarantee that progenitor list is exclusive
ProgenitorData *CrossMatch(Options &opt, const Int_t nbodies, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, long unsigned *&pglist, long unsigned *&noffset, long unsigned *&pfof2, int istepval){
    long int i,j;
    Int_t numshared;
    Double_t merit;
    int nthreads=1,tid;
#ifdef USEOPENMP
#pragma omp parallel 
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    ProgenitorData *p1=new ProgenitorData[nhalos1];
    int *sharelist;
    PriorityQueue *pq;

    if (nhalos2>0){
#ifdef USEOPENMP
    sharelist=new int[nhalos2*nthreads];
    for (i=0;i<nhalos2*nthreads;i++)sharelist[i]=0;
#pragma omp parallel default(shared) \
private(i,j,tid,pq,numshared,merit)
{
#pragma omp for schedule(dynamic)
    for (i=0;i<nhalos1;i++){
        tid=omp_get_thread_num();
        for (j=0;j<h1[i].NumberofParticles;j++){
            if (pfof2[pglist[noffset[i]+j]]>0) sharelist[tid*nhalos2+pfof2[pglist[noffset[i]+j]]-1]+=1;
        }
        numshared=0;
        //determine # of cross-correlations, but if it is not significant relative to Poisson noise, then
        for (j=0;j<nhalos2;j++)
            if(sharelist[tid*nhalos2+j]>opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)){numshared++;}
            else sharelist[tid*nhalos2+j]=0;
        pq=new PriorityQueue(numshared);
        for (j=0;j<nhalos2;j++)if(sharelist[tid*nhalos2+j]>0){
            if (opt.matchtype==NsharedN1N2)
                merit=(Double_t)sharelist[tid*nhalos2+j]*(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            else if (opt.matchtype==NsharedN1)
                merit=(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles;
            else if (opt.matchtype==Nshared)
                merit=(Double_t)sharelist[tid*nhalos2+j];
            else if (opt.matchtype==Nsharedcombo)
                merit=(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[tid*nhalos2+j]*(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            pq->Push(j,merit);
            sharelist[tid*nhalos2+j]=0;
        }
        p1[i].NumberofProgenitors=numshared;
        if (numshared>0) {
            p1[i].ProgenitorList=new long unsigned[numshared];
            p1[i].Merit=new Double_t[numshared];
            p1[i].nsharedfrac=new Double_t[numshared];
            for (j=0;j<numshared;j++){
                //p1[i].ProgenitorList[j]=h2[pq->TopQueue()].haloID;
                p1[i].ProgenitorList[j]=pq->TopQueue();
                p1[i].Merit[j]=pq->TopPriority();
                pq->Pop();
            }
        }
        else {p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;}
        delete pq;
    }
}
    delete[] sharelist;
#else
    sharelist=new int[nhalos2];
    //set max number of progenitors for a single halo
    for (i=0;i<nhalos2;i++)sharelist[i]=0;
    for (i=0;i<nhalos1;i++){
        for (j=0;j<h1[i].NumberofParticles;j++){
            if (pfof2[pglist[noffset[i]+j]]>0) sharelist[pfof2[pglist[noffset[i]+j]]-1]+=1;
        }
        numshared=0;
        //determine # of cross-correlations, but if it is not significant relative to Poisson noise, then
        for (j=0;j<nhalos2;j++)
            if(sharelist[j]>opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)){numshared++;}
            else sharelist[j]=0;
        pq=new PriorityQueue(numshared);
        for (j=0;j<nhalos2;j++)if(sharelist[j]>0){
            if (opt.matchtype==NsharedN1N2)
                merit=(Double_t)sharelist[j]*(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            else if (opt.matchtype==NsharedN1)
                merit=(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles;
            else if (opt.matchtype==Nshared)
                merit=(Double_t)sharelist[j];
            else if (opt.matchtype==Nsharedcombo)
                merit=(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[j]*(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            pq->Push(j,merit);
            sharelist[j]=0;
        }
        p1[i].NumberofProgenitors=numshared;
        if (numshared>0) {
            p1[i].ProgenitorList=new long unsigned[numshared];
            p1[i].Merit=new Double_t[numshared];
            p1[i].nsharedfrac=new Double_t[numshared];
            for (j=0;j<numshared;j++){
                //p1[i].ProgenitorList[j]=h2[pq->TopQueue()].haloID;
                p1[i].ProgenitorList[j]=pq->TopQueue();
                p1[i].Merit[j]=pq->TopPriority();
                pq->Pop();
            }
        }
        else {p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;}
        delete pq;
    }
    delete[] sharelist;
#endif
    }
    else {
    for (i=0;i<nhalos1;i++){
        p1[i].NumberofProgenitors=0;
        p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;
    }
    }
    //adjust number of steps looked back when referencing progenitors
    if (istepval>1) for (i=0;i<nhalos1;i++) p1[i].istep=istepval;

    return p1;
}

///effectively the same code as \ref CrossMatch but allows for the possibility of matching descendant/child nodes using a different merit function
///here the lists are different as input order is reverse of that used in \ref CrossMatch
DescendantData *CrossMatchDescendant(Options &opt, const Int_t nbodies, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, long unsigned *&pglist, long unsigned *&noffset, long unsigned *&pfof2, int istepval){
    long int i,j;
    Int_t numshared;
    Double_t merit;
    Double_t start,end;
    int nthreads=1,tid;
#ifdef USEOPENMP
#pragma omp parallel 
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    DescendantData *d1=new DescendantData[nhalos1];
    int *sharelist;
    PriorityQueue *pq;
    if (nhalos2>0){
#ifdef USEOPENMP
    sharelist=new int[nhalos2*nthreads];
    for (i=0;i<nhalos2*nthreads;i++)sharelist[i]=0;
#pragma omp parallel default(shared) \
private(i,j,tid,pq,numshared,merit)
{
#pragma omp for schedule(dynamic)
    for (i=0;i<nhalos1;i++){
        tid=omp_get_thread_num();
        for (j=0;j<h1[i].NumberofParticles;j++){
            if (pfof2[pglist[noffset[i]+j]]>0) sharelist[tid*nhalos2+pfof2[pglist[noffset[i]+j]]-1]+=1;
        }
        numshared=0;
        //determine # of cross-correlations, but if it is not significant relative to Poisson noise, then
        for (j=0;j<nhalos2;j++)
            if(sharelist[tid*nhalos2+j]>opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)){numshared++;}
            else sharelist[tid*nhalos2+j]=0;
        pq=new PriorityQueue(numshared);
        for (j=0;j<nhalos2;j++)if(sharelist[tid*nhalos2+j]>0){
            if (opt.matchtype==NsharedN1N2)
                merit=(Double_t)sharelist[tid*nhalos2+j]*(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            else if (opt.matchtype==NsharedN1)
                merit=(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles;
            else if (opt.matchtype==Nshared)
                merit=(Double_t)sharelist[tid*nhalos2+j];
            else if (opt.matchtype==Nsharedcombo)
                merit=(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[tid*nhalos2+j]*(Double_t)sharelist[tid*nhalos2+j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            pq->Push(j,merit);
            sharelist[tid*nhalos2+j]=0;
        }
        d1[i].NumberofDescendants=numshared;
        if (numshared>0) {
            d1[i].DescendantList=new long unsigned[numshared];
            d1[i].Merit=new Double_t[numshared];
            d1[i].nsharedfrac=new Double_t[numshared];
            for (j=0;j<numshared;j++){
                //d1[i].DescendantList[j]=h2[pq->TopQueue()].haloID;
                d1[i].DescendantList[j]=pq->TopQueue();
                d1[i].Merit[j]=pq->TopPriority();
                pq->Pop();
            }
        }
        else {d1[i].DescendantList=NULL;d1[i].Merit=NULL;}
        delete pq;
    }
}
    delete[] sharelist;
#else
    sharelist=new int[nhalos2];
    //set max number of progenitors for a single halo
    for (i=0;i<nhalos2;i++)sharelist[i]=0;
    for (i=0;i<nhalos1;i++){
        for (j=0;j<h1[i].NumberofParticles;j++){
            if (pfof2[pglist[noffset[i]+j]]>0) sharelist[pfof2[pglist[noffset[i]+j]]-1]+=1;
            //if (pfof2[pglist1[i][j]]>0) sharelist[pfof2[pglist1[i][j]]-1]+=1;
        }
        numshared=0;
        //determine # of cross-correlations, but if it is not significant relative to Poisson noise, then
        for (j=0;j<nhalos2;j++)
            if(sharelist[j]>opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)){numshared++;}
            else sharelist[j]=0;
        pq=new PriorityQueue(numshared);
        for (j=0;j<nhalos2;j++)if(sharelist[j]>0){
            if (opt.matchtype==NsharedN1N2)
                merit=(Double_t)sharelist[j]*(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            else if (opt.matchtype==NsharedN1)
                merit=(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles;
            else if (opt.matchtype==Nshared)
                merit=(Double_t)sharelist[j];
            else if (opt.matchtype==Nsharedcombo)
                merit=(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[j]*(Double_t)sharelist[j]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
            pq->Push(j,merit);
            sharelist[j]=0;
        }
        d1[i].NumberofDescendants=numshared;
        if (numshared>0) {
            d1[i].DescendantList=new long unsigned[numshared];
            d1[i].Merit=new Double_t[numshared];
            d1[i].nsharedfrac=new Double_t[numshared];
            for (j=0;j<numshared;j++){
                //d1[i].DescendantList[j]=h2[pq->TopQueue()].haloID;
                d1[i].DescendantList[j]=pq->TopQueue();
                d1[i].Merit[j]=pq->TopPriority();
                pq->Pop();
            }
        }
        else {d1[i].DescendantList=NULL;d1[i].Merit=NULL;}
        delete pq;
    }
    delete[] sharelist;
#endif
    }
    else {
    for (i=0;i<nhalos1;i++){
        d1[i].NumberofDescendants=0;
        d1[i].DescendantList=NULL;d1[i].Merit=NULL;
    }
    }
    //adjust number of steps looked forward when referencing descendants
    if (istepval>1) for (i=0;i<nhalos1;i++) d1[i].istep=istepval;
    return d1;
}

//ensure cross match is exclusive
void CleanCrossMatch(const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, ProgenitorData *&p1)
{
    Int_t i,j,k;
    int nthreads=1,tid;
#ifdef USEOPENMP
#pragma omp parallel 
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    int *nh2index=new int[nhalos2];
    int *nh2nummatches=new int[nhalos2];
    //int *count=new int[nhalos2];
    //int **nh2matches=new int*[nhalos2];
    Double_t *merit=new Double_t[nhalos2];
    //store initial matches
    for (i=0;i<nhalos2;i++) {nh2index[i]=merit[i]=-1;nh2nummatches[i]=0;}//count[i]=0;}
    for (i=0;i<nhalos1;i++){
        for(j=0;j<p1[i].NumberofProgenitors;j++) {
            nh2nummatches[p1[i].ProgenitorList[j]]++;
            if(p1[i].Merit[j]>merit[p1[i].ProgenitorList[j]]){
                merit[p1[i].ProgenitorList[j]]=p1[i].Merit[j];
                nh2index[p1[i].ProgenitorList[j]]=i;
            }
        }
    }

    //once stored adjust list ONLY for progenitors with more than 1 descendent
    //for (i=0;i<nhalos2;i++) if (nh2nummatches[i]>=2) {
    //    nh2matches[i]=new int[nh2nummatches[i]];
    //}
    //go through list and if Progenitor's best match is NOT current halo, adjust Progenitor list

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic)
#endif
    for (i=0;i<nhalos1;i++){
        for(j=0;j<p1[i].NumberofProgenitors;j++) {
            if (nh2nummatches[p1[i].ProgenitorList[j]]>=2&&i!=nh2index[p1[i].ProgenitorList[j]]){
                for (k=j;k<p1[i].NumberofProgenitors-1;k++) {
                    p1[i].ProgenitorList[k]=p1[i].ProgenitorList[k+1];
                    p1[i].Merit[k]=p1[i].Merit[k+1];
                }
                j--;p1[i].NumberofProgenitors--;
            }
        }
    }
#ifdef USEOPENMP
}
#endif

    //adjust data to store haloIDS
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic)
#endif
    for (i=0;i<nhalos1;i++){
        for(j=0;j<p1[i].NumberofProgenitors;j++) {
            //store nshared
            p1[i].nsharedfrac[j]=sqrt(p1[i].Merit[j]*((double)h2[p1[i].ProgenitorList[j]].NumberofParticles*(double)h1[i].NumberofParticles))/(double)h1[i].NumberofParticles;
            p1[i].ProgenitorList[j]=h2[p1[i].ProgenitorList[j]].haloID;
        }
    }
#ifdef USEOPENMP
}
#endif
    delete[] nh2index;
    delete[] nh2nummatches;
    delete[] merit;
}

void UpdateRefProgenitors(const Int_t numhalos, ProgenitorData *&pref, ProgenitorData *&ptemp)
{
    for (Int_t i=0;i<numhalos;i++) 
        if (pref[i].NumberofProgenitors==0 && ptemp[i].NumberofProgenitors>0) pref[i]=ptemp[i];
}

void CleanCrossMatchDescendant(const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, DescendantData *&p1)
{
    Int_t i,j,k;
    int nthreads=1,tid;
#ifdef USEOPENMP
#pragma omp parallel 
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    int *nh2index=new int[nhalos2];
    int *nh2nummatches=new int[nhalos2];
    //int *count=new int[nhalos2];
    //int **nh2matches=new int*[nhalos2];
    Double_t *merit=new Double_t[nhalos2];
    //store initial matches
    for (i=0;i<nhalos2;i++) {nh2index[i]=merit[i]=-1;nh2nummatches[i]=0;}//count[i]=0;}
    for (i=0;i<nhalos1;i++){
        for(j=0;j<p1[i].NumberofDescendants;j++) {
            nh2nummatches[p1[i].DescendantList[j]]++;
            if(p1[i].Merit[j]>merit[p1[i].DescendantList[j]]){
                merit[p1[i].DescendantList[j]]=p1[i].Merit[j];
                nh2index[p1[i].DescendantList[j]]=i;
            }
        }
    }

    //once stored adjust list ONLY for progenitors with more than 1 descendent
    //for (i=0;i<nhalos2;i++) if (nh2nummatches[i]>=2) {
    //    nh2matches[i]=new int[nh2nummatches[i]];
    //}
    //go through list and if Progenitor's best match is NOT current halo, adjust Progenitor list

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic)
#endif
    for (i=0;i<nhalos1;i++){
        for(j=0;j<p1[i].NumberofDescendants;j++) {
            if (nh2nummatches[p1[i].DescendantList[j]]>=2&&i!=nh2index[p1[i].DescendantList[j]]){
                for (k=j;k<p1[i].NumberofDescendants-1;k++) {
                    p1[i].DescendantList[k]=p1[i].DescendantList[k+1];
                    p1[i].Merit[k]=p1[i].Merit[k+1];
                }
                j--;p1[i].NumberofDescendants--;
            }
        }
    }
#ifdef USEOPENMP
}
#endif

    //adjust data to store haloIDS
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic)
#endif
    for (i=0;i<nhalos1;i++){
        for(j=0;j<p1[i].NumberofDescendants;j++) {
            //store nshared
            p1[i].nsharedfrac[j]=sqrt(p1[i].Merit[j]*((double)h2[p1[i].DescendantList[j]].NumberofParticles*(double)h1[i].NumberofParticles))/(double)h1[i].NumberofParticles;
            p1[i].DescendantList[j]=h2[p1[i].DescendantList[j]].haloID;
        }
    }
#ifdef USEOPENMP
}
#endif
    delete[] nh2index;
    delete[] nh2nummatches;
    delete[] merit;
}

void UpdateRefDescendants(const Int_t numhalos, DescendantData *&dref, DescendantData *&dtemp)
{
    for (Int_t i=0;i<numhalos;i++) 
        if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) dref[i]=dtemp[i];
}

//@}

/// \name if particle ids need to be mapped to indices
//@{
void MapPIDStoIndex(Options &opt, HaloTreeData *&pht) {
    cout<<"Mapping PIDS to index "<<endl;
    Int_t i,j,k;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<opt.numsnapshots;i++) {
            for (j=0;j<pht[i].numhalos;j++) 
                for (k=0;k<pht[i].Halo[j].NumberofParticles;k++) 
                    opt.mappingfunc(pht[i].Halo[j].ParticleID[k]);
    }
#ifdef USEOPENMP
}
#endif
}

void simplemap(long unsigned int &i) {}

//@}

/// \name Simple group id based array building
//@{

///build group size array
Int_t *BuildNumInGroup(const Int_t nbodies, const Int_t numgroups, Int_t *pfof){
    Int_t *numingroup=new Int_t[numgroups+1];
    for (Int_t i=0;i<=numgroups;i++) numingroup[i]=0;
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) numingroup[pfof[i]]++;
    return numingroup;
}
///build the group particle index list (assumes particles are in file Index order)
Int_t **BuildPGList(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t *pfof){
    Int_t **pglist=new Int_t*[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++) {pglist[i]=new Int_t[numingroup[i]];numingroup[i]=0;}
    //for (Int_t i=0;i<nbodies;i++) if (ids[i]!=-1) {pglist[pfof[ids[i]]][numingroup[pfof[ids[i]]]++]=ids[i];}
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) pglist[pfof[i]][numingroup[pfof[i]]++]=i;
    return pglist;
}
///build the Head array which points to the head of the group a particle belongs to
Int_t *BuildHeadArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_t *Head=new Int_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) Head[i]=i;
    for (Int_t i=1;i<=numgroups;i++)
        for (Int_t j=1;j<numingroup[i];j++) Head[pglist[i][j]]=Head[pglist[i][0]];
    return Head;
}
///build the Next array which points to the next particle in the group
Int_t *BuildNextArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_t *Next=new Int_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) Next[i]=-1;
    for (Int_t i=1;i<=numgroups;i++)
        for (Int_t j=0;j<numingroup[i]-1;j++) Next[pglist[i][j]]=pglist[i][j+1];
    return Next;
}
///build the Len array which stores the length of the group a particle belongs to 
Int_t *BuildLenArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_t *Len=new Int_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) Len[i]=0;
    for (Int_t i=1;i<=numgroups;i++)
        for (Int_t j=0;j<numingroup[i];j++) Len[pglist[i][j]]=numingroup[i];
    return Len;
}
///build the GroupTail array which stores the Tail of a group
Int_t *BuildGroupTailArray(const Int_t nbodies, const Int_t numgroups, Int_t *numingroup, Int_t **pglist){
    Int_t *GTail=new Int_t[numgroups+1];
    for (Int_t i=1;i<=numgroups;i++)
        GTail[i]=pglist[i][numingroup[i]-1];
    return GTail;
}
//@}
