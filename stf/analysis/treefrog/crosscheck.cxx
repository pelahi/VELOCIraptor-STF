
/*! \file crosscheck.cxx
 *  \brief this file contains routines for determining and calculating cross matches
 */

#include "TreeFrog.h"


/// \name cross matching routines
//@{

/// Determine initial progenitor list based on just merit
/// does not guarantee that progenitor list is exclusive
ProgenitorData *CrossMatch(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, unsigned int*&pfof2, int &ilistupdated, int istepval, ProgenitorData *refprogen)
{
    long int i,j,k,n,index;
    Int_t numshared;
    Double_t merit;
    int nthreads=1,tid,chunksize;
    long unsigned offset;
    //temp variable to store the openmp reduction value for ilistupdated
    int newilistupdated;
#ifdef USEOPENMP
#pragma omp parallel 
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
    //chunksize should be small as typically halo distribution is a few very large haloes followed by many smaller ones
    //to ensure reasonably even splitting of work, have small chunks
    chunksize=min((long unsigned)(nhalos1/nthreads+1),OMPCHUNKSIZE);
#endif
    ProgenitorData *p1=new ProgenitorData[nhalos1];
    int *sharelist,*halolist;
    PriorityQueue *pq,*pq2;
    int np1,np2;
    long unsigned num_noprogen, ntotitems;
    long unsigned *needprogenlist;
    //init that list is updated if no reference list is provided
    if (refprogen==NULL) ilistupdated=1;
    //otherwise assume list is not updated
    else ilistupdated=0;
    newilistupdated=ilistupdated;

    //if there are halos to link
    if (nhalos2>0){
    //if a reference list is not provided, then build for every halo
    if (refprogen==NULL) {
    ntotitems=nhalos2*(long unsigned)nthreads;
    sharelist=new int[ntotitems];
    halolist=new int[ntotitems];
    //to store haloes that share links and the halo index of those shared haloes
    for (i=0;i<ntotitems;i++)sharelist[i]=0;
#ifdef USEOPENMP
#pragma omp parallel for schedule(dynamic,chunksize) \
default(shared) \
private(j,k,tid,pq,numshared,merit,index,offset,np1,np2,pq2)
#endif
    for (i=0;i<nhalos1;i++){
#ifdef USEOPENMP
        //initialize variables
        tid=omp_get_thread_num();
        offset=((long int)tid)*nhalos2;
#else
        offset=0;
#endif
        numshared=0;
        //go through halos particle list to see if these particles belong to another halo
        //at a different time/in a different catalog
        for (j=0;j<h1[i].NumberofParticles;j++){
            index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
            if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
            //if first time halo has been added, update halolist and increase numshared
            if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
        }
        //now proccess numshared list to remove insignificant connections
        for (k=0;k<numshared;k++) {
            j=halolist[offset+k];
            index=offset+j;
            //if sharelist not enough, then move to end and decrease numshared
            if(sharelist[index]<opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)) {
                halolist[offset+k]=halolist[offset+numshared-1];
                sharelist[index]=0;
                k--;
                numshared--;
            }
        }
        //now calculate merits and sort according to best merit.
        if (numshared>0) {
            //store all viable matches
            pq=new PriorityQueue(numshared);
            for (k=0;k<numshared;k++) {
                j=halolist[offset+k];
                index=offset+j;
                if (opt.matchtype==NsharedN1N2)
                    merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                else if (opt.matchtype==NsharedN1)
                    merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles;
                else if (opt.matchtype==Nshared)
                    merit=(Double_t)sharelist[index];
                else if (opt.matchtype==Nsharedcombo)
                    merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                pq->Push(j,merit);
                sharelist[index]=0;
            }
            p1[i].NumberofProgenitors=numshared;
            p1[i].ProgenitorList=new long unsigned[numshared];
            p1[i].Merit=new float[numshared];
            p1[i].nsharedfrac=new float[numshared];
            for (j=0;j<numshared;j++){
                //p1[i].ProgenitorList[j]=h2[pq->TopQueue()].haloID;
                p1[i].ProgenitorList[j]=pq->TopQueue();
                p1[i].Merit[j]=pq->TopPriority();
                pq->Pop();
            }
            delete pq;
        }
        //if the number shared is zero, do nothing
        else {p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;}

        //if weighted merit function is to be calculated then use the most bound fraction of particles to construct share list
        if (opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
            //calculate number of particles 
            np1=(h1[i].NumberofParticles*opt.particle_frac);
            //if halo has enough particle for meaningful most bound particles check, proceed to find merits. Similar to normal merit 
            //calculations but have now np2 for calculating merits
            if (np1>=opt.min_numpart) { 
                numshared=0;
                for (j=0;j<np1;j++){
                    index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
                    if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
                    if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
                }
                for (k=0;k<numshared;k++) {
                    j=halolist[offset+k];
                    index=offset+j;
                    np2=h2[j].NumberofParticles*opt.particle_frac;
                    if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                    else if (np2<opt.min_numpart) np2=opt.min_numpart;
                    //if sharelist not enough, then move to end and decrease numshared
                    if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2)){
                        halolist[offset+k]=halolist[offset+numshared-1];
                        sharelist[index]=0;
                        k--;
                        numshared--;
                    }
                }
                if (numshared>0) {
                    //store all viable matches and old ones too
                    pq2=new PriorityQueue(numshared+p1[i].NumberofProgenitors);
                    for (k=0;k<numshared;k++) {
                        j=halolist[offset+k];
                        index=offset+j;
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
                        if (opt.matchtype==NsharedN1N2)
                            merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        else if (opt.matchtype==NsharedN1)
                            merit=(Double_t)sharelist[index]/(Double_t)np1;
                        else if (opt.matchtype==Nshared)
                            merit=(Double_t)sharelist[index];
                        else if (opt.matchtype==Nsharedcombo)
                            merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        pq2->Push(j,merit);
                        sharelist[index]=0;
                    }
                    //really only care about adjust the best match found if any
                    for (j=0;j<p1[i].NumberofProgenitors;j++) pq2->Push(p1[i].ProgenitorList[j],p1[i].Merit[j]);
                    //if most bound particle matching is not the same as overall matching, lets rearrange match
                    if (pq2->TopQueue()!=p1[i].ProgenitorList[0]) {
                        for (j=0;j<p1[i].NumberofProgenitors;j++) if (p1[i].ProgenitorList[j]==pq2->TopQueue()) {
                            p1[i].ProgenitorList[j]=p1[i].ProgenitorList[0];
                            p1[i].Merit[j]=p1[i].Merit[0];
                            p1[i].ProgenitorList[0]=pq2->TopQueue();
                            p1[i].Merit[0]=pq2->TopPriority();
                            break;
                        }
                        delete pq2;
                    }
                }
                //end of numshared>0 check
            }
            //end of halo has enough particles for weighted (bound fraction) merit adjustment
        }
        //end of weighted (based on most bound particles) merit 
    }
    delete[] sharelist;
    delete[] halolist;
    }
    //if a reference list is provided then only link halos with no current link
    else {
    //check to see if haloes need to search at all
    num_noprogen=0;
    for (i=0;i<nhalos1;i++){
        p1[i].NumberofProgenitors=0;p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;
        if (refprogen[i].NumberofProgenitors==0) {
            num_noprogen+=1;
        }
        else {
            if (refprogen[i].Merit[0]<opt.meritlimit) num_noprogen+=1;
        }
    }
    //only allocate memory and process list if there are any haloes needing to be searched
    if (num_noprogen>0) {
        needprogenlist=new long unsigned[num_noprogen];
        num_noprogen=0;
        for (i=0;i<nhalos1;i++){
            if (refprogen[i].NumberofProgenitors==0){
                needprogenlist[num_noprogen++]=i;
            }
            else {
                if (refprogen[i].Merit[0]<opt.meritlimit) needprogenlist[num_noprogen++]=i;
            }
        }

        ntotitems=nhalos2*(long unsigned)nthreads;
        sharelist=new int[ntotitems];
        halolist=new int[ntotitems];
        for (i=0;i<ntotitems;i++)sharelist[i]=0;

#ifdef USEOPENMP
#pragma omp parallel for schedule(dynamic,chunksize) reduction(+:newilistupdated) \
default(shared) \
private(i,j,n,tid,pq,numshared,merit,index,offset,np1,np2,pq2)
#endif
        for (k=0;k<num_noprogen;k++){
#ifdef USEOPENMP
            tid=omp_get_thread_num();
            offset=((long int)tid)*nhalos2;
#else
            offset=0;
#endif
            i=needprogenlist[k];
            numshared=0;
            for (j=0;j<h1[i].NumberofParticles;j++){
                index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
                if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
                //if first time halo has been added, update halolist and increase numshared
                if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
            }
            //now proccess numshared list to remove insignificant connections
            for (n=0;n<numshared;n++) {
                j=halolist[offset+n];
                index=offset+j;
                //if sharelist not enough, then move to end and decrease numshared
                if(sharelist[index]<opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)) {
                    halolist[offset+n]=halolist[offset+numshared-1];
                    sharelist[index]=0;
                    n--;
                    numshared--;
                }
            }
            if (numshared>0) {
                pq=new PriorityQueue(numshared);
                for (n=0;n<numshared;n++){
                    j=halolist[offset+n];
                    index=offset+j;
                    if (opt.matchtype==NsharedN1N2)
                        merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                    else if (opt.matchtype==NsharedN1)
                        merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles;
                    else if (opt.matchtype==Nshared)
                        merit=(Double_t)sharelist[index];
                    else if (opt.matchtype==Nsharedcombo)
                        merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                    pq->Push(j,merit);
                    sharelist[index]=0;
                }
                //if at this point when looking for progenitors not present in the reference list
                //can check to see if numshared>0 and if so, then reference list will have to be updated
                newilistupdated+=(numshared>0);

                p1[i].NumberofProgenitors=numshared;
                p1[i].ProgenitorList=new long unsigned[numshared];
                p1[i].Merit=new float[numshared];
                p1[i].nsharedfrac=new float[numshared];
                for (j=0;j<numshared;j++){
                    //p1[i].ProgenitorList[j]=h2[pq->TopQueue()].haloID;
                    p1[i].ProgenitorList[j]=pq->TopQueue();
                    p1[i].Merit[j]=pq->TopPriority();
                    pq->Pop();
                }
                delete pq;
            }
            else {p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;}
            //if weighted merit function is to be calculated then use the most bound fraction of particles to construct share list
            if (opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
                np1=(h1[i].NumberofParticles*opt.particle_frac);
                if (np1>=opt.min_numpart) { 
                    numshared=0;
                    for (j=0;j<np1;j++){
                        index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
                        if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
                        if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
                    }
                    for (n=0;n<numshared;n++) {
                        j=halolist[offset+n];
                        index=offset+j;
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
                        //if sharelist not enough, then move to end and decrease numshared
                        if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2)){
                            halolist[offset+n]=halolist[offset+numshared-1];
                            sharelist[index]=0;
                            n--;
                            numshared--;
                        }
                    }
                    if (numshared>0) {
                        //store all viable matches and old ones too
                        pq2=new PriorityQueue(numshared+p1[i].NumberofProgenitors);
                        for (n=0;n<numshared;n++) {
                            j=halolist[offset+n];
                            index=offset+j;
                            np2=h2[j].NumberofParticles*opt.particle_frac;
                            if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                            else if (np2<opt.min_numpart) np2=opt.min_numpart;
                            if (opt.matchtype==NsharedN1N2)
                                merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                            else if (opt.matchtype==NsharedN1)
                                merit=(Double_t)sharelist[index]/(Double_t)np1;
                            else if (opt.matchtype==Nshared)
                                merit=(Double_t)sharelist[index];
                            else if (opt.matchtype==Nsharedcombo)
                                merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                            pq2->Push(j,merit);
                            sharelist[index]=0;
                        }
                        //really only care about adjust the best match found if any
                        for (j=0;j<p1[i].NumberofProgenitors;j++) pq2->Push(p1[i].ProgenitorList[j],p1[i].Merit[j]);
                        //if most bound particle matching is not the same as overall matching, lets rearrange match
                        if (pq2->TopQueue()!=p1[i].ProgenitorList[0]) {
                            for (j=0;j<p1[i].NumberofProgenitors;j++) if (p1[i].ProgenitorList[j]==pq2->TopQueue()) {
                                p1[i].ProgenitorList[j]=p1[i].ProgenitorList[0];
                                p1[i].Merit[j]=p1[i].Merit[0];
                                p1[i].ProgenitorList[0]=pq2->TopQueue();
                                p1[i].Merit[0]=pq2->TopPriority();
                                break;
                            }
                            delete pq2;
                        }
                    }
                }
            }
            //end of weighted merit
        }
        delete[] sharelist;
        delete[] halolist;
        delete[] needprogenlist;
    }
    //end of if statement for progenitors more than one snap ago if any current haloes need to be searched
    ilistupdated=newilistupdated;
    }
    //end of progenitors more than one snap ago
    }
    //if no halos then there are no links
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
DescendantData *CrossMatchDescendant(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, unsigned int *&pfof2, int &ilistupdated, int istepval, DescendantData *refdescen)
{
    long int i,j,k,n,index;
    Int_t numshared;
    Double_t merit;
    Double_t start,end;
    int nthreads=1,tid;
    long unsigned offset;
    //temp variable to store the openmp reduction value for ilistupdated
    int newilistupdated;
    int chunksize;
#ifdef USEOPENMP
#pragma omp parallel 
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
    //chunksize should be small as typically halo distribution is a few very large haloes followed by many smaller ones
    //to ensure reasonably even splitting of work, have small chunks
    chunksize=min((long unsigned)(nhalos1/nthreads+1),OMPCHUNKSIZE);
#endif
    DescendantData *d1=new DescendantData[nhalos1];
    int *sharelist, *halolist;
    PriorityQueue *pq,*pq2;
    int np1,np2;
    long unsigned num_nodescen, ntotitems;
    long unsigned *needdescenlist;
    if (refdescen==NULL) ilistupdated=1;
    else ilistupdated=0;
    newilistupdated=ilistupdated;

    if (nhalos2>0){

    if (refdescen==NULL) {
    ntotitems=nhalos2*(long unsigned)nthreads;
    sharelist=new int[ntotitems];
    halolist=new int[ntotitems];
    //to store haloes that share links and the halo index of those shared haloes
    for (i=0;i<ntotitems;i++)sharelist[i]=0;
#ifdef USEOPENMP
#pragma omp parallel for schedule(dynamic,chunksize) \
default(shared) \
private(i,j,k,tid,pq,numshared,merit,index,offset,np1,np2,pq2)
#endif
    for (i=0;i<nhalos1;i++){
#ifdef USEOPENMP
        //initialize variables
        tid=omp_get_thread_num();
        offset=((long int)tid)*nhalos2;
#else
        offset=0;
#endif
        numshared=0;
        //go through halos particle list to see if these particles belong to another halo
        //at a different time/in a different catalog
        for (j=0;j<h1[i].NumberofParticles;j++){
            index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
            if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
            //if first time halo has been added, update halolist and increase numshared
            if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
        }
        //now proccess numshared list to remove insignificant connections
        for (k=0;k<numshared;k++) {
            j=halolist[offset+k];
            index=offset+j;
            //if sharelist not enough, then move to end and decrease numshared
            if(sharelist[index]<opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)) {
                halolist[offset+k]=halolist[offset+numshared-1];
                sharelist[index]=0;
                k--;
                numshared--;
            }
        }
        //now calculate merits and sort according to best merit.
        if (numshared>0) {
            //store all viable matches
            pq=new PriorityQueue(numshared);
            for (k=0;k<numshared;k++) {
                j=halolist[offset+k];
                index=offset+j;
                if (opt.matchtype==NsharedN1N2)
                    merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                else if (opt.matchtype==NsharedN1)
                    merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles;
                else if (opt.matchtype==Nshared)
                    merit=(Double_t)sharelist[index];
                else if (opt.matchtype==Nsharedcombo)
                    merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                pq->Push(j,merit);
                sharelist[index]=0;
            }
            d1[i].NumberofDescendants=numshared;
            d1[i].DescendantList=new long unsigned[numshared];
            d1[i].Merit=new float[numshared];
            d1[i].nsharedfrac=new float[numshared];
            for (j=0;j<numshared;j++){
                //d1[i].DescendantList[j]=h2[pq->TopQueue()].haloID;
                d1[i].DescendantList[j]=pq->TopQueue();
                d1[i].Merit[j]=pq->TopPriority();
                pq->Pop();
            }
            delete pq;
        }
        //if the number shared is zero, do nothing
        else {d1[i].DescendantList=NULL;d1[i].Merit=NULL;}

        //if weighted merit function is to be calculated then use the most bound fraction of particles to construct share list
        if (opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
            //calculate number of particles 
            np1=(h1[i].NumberofParticles*opt.particle_frac);
            //if halo has enough particle for meaningful most bound particles check, proceed to find merits. Similar to normal merit 
            //calculations but have now np2 for calculating merits
            if (np1>=opt.min_numpart) { 
                numshared=0;
                for (j=0;j<np1;j++){
                    index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
                    if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
                    if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
                }
                for (k=0;k<numshared;k++) {
                    j=halolist[offset+k];
                    index=offset+j;
                    np2=h2[j].NumberofParticles*opt.particle_frac;
                    if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                    else if (np2<opt.min_numpart) np2=opt.min_numpart;
                    //if sharelist not enough, then move to end and decrease numshared
                    if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2)){
                        halolist[offset+k]=halolist[offset+numshared-1];
                        sharelist[index]=0;
                        k--;
                        numshared--;
                    }
                }
                if (numshared>0) {
                    //store all viable matches and old ones too
                    pq2=new PriorityQueue(numshared+d1[i].NumberofDescendants);
                    for (k=0;k<numshared;k++) {
                        j=halolist[offset+k];
                        index=offset+j;
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
                        if (opt.matchtype==NsharedN1N2)
                            merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        else if (opt.matchtype==NsharedN1)
                            merit=(Double_t)sharelist[index]/(Double_t)np1;
                        else if (opt.matchtype==Nshared)
                            merit=(Double_t)sharelist[index];
                        else if (opt.matchtype==Nsharedcombo)
                            merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        pq2->Push(j,merit);
                        sharelist[index]=0;
                    }
                    //really only care about adjust the best match found if any
                    for (j=0;j<d1[i].NumberofDescendants;j++) pq2->Push(d1[i].DescendantList[j],d1[i].Merit[j]);
                    //if most bound particle matching is not the same as overall matching, lets rearrange match
                    if (pq2->TopQueue()!=d1[i].DescendantList[0]) {
                        for (j=0;j<d1[i].NumberofDescendants;j++) if (d1[i].DescendantList[j]==pq2->TopQueue()) {
                            d1[i].DescendantList[j]=d1[i].DescendantList[0];
                            d1[i].Merit[j]=d1[i].Merit[0];
                            d1[i].DescendantList[0]=pq2->TopQueue();
                            d1[i].Merit[0]=pq2->TopPriority();
                            break;
                        }
                        delete pq2;
                    }
                }
                //end of numshared>0 check
            }
            //end of halo has enough particles for weighted (bound fraction) merit adjustment
        }
        //end of weighted (based on most bound particles) merit 
    }
        delete[] sharelist;
        delete[] halolist;
    }
    else {
    //check to see if haloes need to search at all
    num_nodescen=0;
    for (i=0;i<nhalos1;i++){
        d1[i].NumberofDescendants=0;d1[i].DescendantList=NULL;d1[i].Merit=NULL;
        if (refdescen[i].NumberofDescendants==0) {
            num_nodescen+=1;
        }
        else {
            if (refdescen[i].Merit[0]<opt.meritlimit) num_nodescen+=1;
        }
    }
    //only allocate memory and process list if there are any haloes needing to be searched
    if (num_nodescen>0) {
        needdescenlist=new long unsigned[num_nodescen];
        num_nodescen=0;
        for (i=0;i<nhalos1;i++){
            if (refdescen[i].NumberofDescendants==0){
                needdescenlist[num_nodescen++]=i;
            }
            else {
                if (refdescen[i].Merit[0]<opt.meritlimit) needdescenlist[num_nodescen++]=i;
            }
        }

        ntotitems=nhalos2*(long unsigned)nthreads;
        sharelist=new int[ntotitems];
        halolist=new int[ntotitems];
        for (i=0;i<ntotitems;i++)sharelist[i]=0;

#ifdef USEOPENMP
#pragma omp parallel for schedule(dynamic,chunksize) reduction(+:newilistupdated) \
default(shared) \
private(i,j,n,tid,pq,numshared,merit,index,offset,np1,np2,pq2)
#endif
        for (k=0;k<num_nodescen;k++){
#ifdef USEOPENMP
            tid=omp_get_thread_num();
            offset=((long int)tid)*nhalos2;
#else
            offset=0;
#endif
            i=needdescenlist[k];
            numshared=0;
            for (j=0;j<h1[i].NumberofParticles;j++){
                index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
                if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
                //if first time halo has been added, update halolist and increase numshared
                if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
            }
            //now proccess numshared list to remove insignificant connections
            for (n=0;n<numshared;n++) {
                j=halolist[offset+n];
                index=offset+j;
                //if sharelist not enough, then move to end and decrease numshared
                if(sharelist[index]<opt.mlsig*sqrt((Double_t)h2[j].NumberofParticles)) {
                    halolist[offset+n]=halolist[offset+numshared-1];
                    sharelist[index]=0;
                    n--;
                    numshared--;
                }
            }
            if (numshared>0) {
                pq=new PriorityQueue(numshared);
                for (n=0;n<numshared;n++){
                    j=halolist[offset+n];
                    index=offset+j;
                    if (opt.matchtype==NsharedN1N2)
                        merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                    else if (opt.matchtype==NsharedN1)
                        merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles;
                    else if (opt.matchtype==Nshared)
                        merit=(Double_t)sharelist[index];
                    else if (opt.matchtype==Nsharedcombo)
                        merit=(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)h1[i].NumberofParticles/(Double_t)h2[j].NumberofParticles;
                    pq->Push(j,merit);
                    sharelist[index]=0;
                }
                //if at this point when looking for progenitors not present in the reference list
                //can check to see if numshared>0 and if so, then reference list will have to be updated
                newilistupdated+=(numshared>0);

                d1[i].NumberofDescendants=numshared;
                d1[i].DescendantList=new long unsigned[numshared];
                d1[i].Merit=new float[numshared];
                d1[i].nsharedfrac=new float[numshared];
                for (j=0;j<numshared;j++){
                    //d1[i].DescendantList[j]=h2[pq->TopQueue()].haloID;
                    d1[i].DescendantList[j]=pq->TopQueue();
                    d1[i].Merit[j]=pq->TopPriority();
                    pq->Pop();
                }
                delete pq;
            }
            else {d1[i].DescendantList=NULL;d1[i].Merit=NULL;}
            //if weighted merit function is to be calculated then use the most bound fraction of particles to construct share list
            if (opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
                np1=(h1[i].NumberofParticles*opt.particle_frac);
                if (np1>=opt.min_numpart) { 
                    numshared=0;
                    for (j=0;j<np1;j++){
                        index=offset+pfof2[h1[i].ParticleID[j]]-(long int)1;
                        if (pfof2[h1[i].ParticleID[j]]>0) sharelist[index]+=1;
                        if (sharelist[index]==1) halolist[offset+numshared++]=pfof2[h1[i].ParticleID[j]]-1;
                    }
                    for (n=0;n<numshared;n++) {
                        j=halolist[offset+n];
                        index=offset+j;
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
                        //if sharelist not enough, then move to end and decrease numshared
                        if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2)){
                            halolist[offset+n]=halolist[offset+numshared-1];
                            sharelist[index]=0;
                            n--;
                            numshared--;
                        }
                    }
                    if (numshared>0) {
                        //store all viable matches and old ones too
                        pq2=new PriorityQueue(numshared+d1[i].NumberofDescendants);
                        for (n=0;n<numshared;n++) {
                            j=halolist[offset+n];
                            index=offset+j;
                            np2=h2[j].NumberofParticles*opt.particle_frac;
                            if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                            else if (np2<opt.min_numpart) np2=opt.min_numpart;
                            if (opt.matchtype==NsharedN1N2)
                                merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                            else if (opt.matchtype==NsharedN1)
                                merit=(Double_t)sharelist[index]/(Double_t)np1;
                            else if (opt.matchtype==Nshared)
                                merit=(Double_t)sharelist[index];
                            else if (opt.matchtype==Nsharedcombo)
                                merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                            pq2->Push(j,merit);
                            sharelist[index]=0;
                        }
                        //really only care about adjust the best match found if any
                        for (j=0;j<d1[i].NumberofDescendants;j++) pq2->Push(d1[i].DescendantList[j],d1[i].Merit[j]);
                        //if most bound particle matching is not the same as overall matching, lets rearrange match
                        if (pq2->TopQueue()!=d1[i].DescendantList[0]) {
                            for (j=0;j<d1[i].NumberofDescendants;j++) if (d1[i].DescendantList[j]==pq2->TopQueue()) {
                                d1[i].DescendantList[j]=d1[i].DescendantList[0];
                                d1[i].Merit[j]=d1[i].Merit[0];
                                d1[i].DescendantList[0]=pq2->TopQueue();
                                d1[i].Merit[0]=pq2->TopPriority();
                                break;
                            }
                            delete pq2;
                        }
                    }
                }
            }
            //end of weighted merit
        }
        delete[] sharelist;
        delete[] halolist;
        delete[] needdescenlist;
    }
    //end of if statement for progenitors more than one snap ago if any current haloes need to be searched
    ilistupdated=newilistupdated;
    }
    //end of progenitors more than one snap ago
    }
    //end of more than zero halos to produce links to
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

///ensure cross match is exclusive
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
#pragma omp for schedule(dynamic,10) nowait
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
#pragma omp for schedule(dynamic,10) nowait
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
///similar to \ref CleanCrossMatch but does the same for descendants
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
#pragma omp for schedule(dynamic,10) nowait
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
#pragma omp for schedule(dynamic,10) nowait
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
///when linking using multiple snapshots, use to update the progenitor list based on a candidate list built using a single snapshot
void UpdateRefProgenitors(const int ilink, const Int_t numhalos, ProgenitorData *&pref, ProgenitorData *&ptemp)
{
    if (ilink==MSLCMISSING) {
        for (Int_t i=0;i<numhalos;i++) 
            if (pref[i].NumberofProgenitors==0 && ptemp[i].NumberofProgenitors>0) pref[i]=ptemp[i];
    }
    else if (ilink==MSLCMERIT) {
        for (Int_t i=0;i<numhalos;i++) {
            if (pref[i].NumberofProgenitors==0 && ptemp[i].NumberofProgenitors>0) pref[i]=ptemp[i];
            else if (pref[i].NumberofProgenitors>0 && ptemp[i].NumberofProgenitors>0 && pref[i].Merit[0]<ptemp[i].Merit[0]) pref[i]=ptemp[i];
        }
    }
}
///similar to \ref UpdateRefProgenitors but for descendants
void UpdateRefDescendants(const int ilink, const Int_t numhalos, DescendantData *&dref, DescendantData *&dtemp)
{
    if (ilink==MSLCMISSING) {
        for (Int_t i=0;i<numhalos;i++) 
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) dref[i]=dtemp[i];
    }
    else if (ilink==MSLCMERIT) {
        for (Int_t i=0;i<numhalos;i++) {
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) dref[i]=dtemp[i];
            else if (dref[i].NumberofDescendants>0 && dtemp[i].NumberofDescendants>0 && dref[i].Merit[0]<dtemp[i].Merit[0]) dref[i]=dtemp[i];
        }
    }
}

///Construction of descendant list using the candidate progenitor list
void BuildProgenitorBasedDescendantList(Int_t itimeprogen, Int_t itimedescen, Int_t nhalos, ProgenitorData *&pprogen, DescendantDataProgenBased **&pprogendescen, int istep)
{
    //if first pass then store descendants looking at all progenitors
    int did;
    for (Int_t k=0;k<nhalos;k++) {
        for (Int_t nprogs=0;nprogs<pprogen[k].NumberofProgenitors;nprogs++) if (pprogen[k].istep==istep) {
            did=pprogen[k].ProgenitorList[nprogs]-1;//make sure halo descendent index set to start at 0
            pprogendescen[itimedescen][did].haloindex.push_back(k);
            pprogendescen[itimedescen][did].halotemporalindex.push_back(itimeprogen);
            pprogendescen[itimedescen][did].Merit.push_back(pprogen[k].Merit[nprogs]);
            pprogendescen[itimedescen][did].deltat.push_back(istep);
#ifdef USEMPI
            pprogendescen[itimedescen][did].MPITask.push_back(ThisTask);
#endif
            pprogendescen[itimedescen][did].NumberofDescendants++;
        }
    }
}

///Cleans the progenitor list to ensure that a given object only has a single descendant using the descendant list built using progenitors. Called if linking across multiple snapshots
void CleanProgenitorsUsingDescendants(Int_t i, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen)
{
    unsigned long itime;
    unsigned long hid;
#ifdef USEMPI
    int itask;
#endif
    int itemp;
    //parse list of descendants for objects with more than one and adjust the progenitor list
    for (Int_t k=0;k<pht[i].numhalos;k++) if (pprogendescen[i][k].NumberofDescendants>1){
        //adjust the progenitor list of all descendents affected, ie: that will have the progenitor removed from the list
        //note that the optimal descendant is always placed at position 0
        pprogendescen[i][k].OptimalTemporalMerit();
#ifdef USEMPI
        //if mpi, it is possible that the best match is not a local halo 
        ///\todo do I need to do anything when the best descendant is not local to mpi domain? I don't think so since only 
        ///matters if I am changing the local list
#endif
        for (Int_t n=1;n<pprogendescen[i][k].NumberofDescendants;n++) {
            itime=pprogendescen[i][k].halotemporalindex[n];
            hid=pprogendescen[i][k].haloindex[n];
#ifdef USEMPI
            //if mpi is invoked, it possible the progenitor that is supposed to have a descendant removed is not on this mpi task 
            itask=pprogendescen[i][k].MPITask[n];
            //only change progenitors of halos on this task
            if (itask==ThisTask) {
#endif
            itemp=0;
            //k+1 is halo indexing versus halo id values
            while (pprogen[itime][hid].ProgenitorList[itemp]!=k+1) {itemp++;}
            for (Int_t j=itemp;j<pprogen[itime][hid].NumberofProgenitors-1;j++) {
                pprogen[itime][hid].ProgenitorList[j]=pprogen[itime][hid].ProgenitorList[j+1];
                pprogen[itime][hid].Merit[j]=pprogen[itime][hid].Merit[j+1];
                if (pprogen[itime][hid].nsharedfrac!=NULL) pprogen[itime][hid].nsharedfrac[j]=pprogen[itime][hid].nsharedfrac[j+1];
            }
            pprogen[itime][hid].NumberofProgenitors--;
#ifdef USEMPI
            }
#endif
        }
    }
}

//@}

/// \name if halo ids need to be adjusted
//@{
void UpdateHaloIDs(Options &opt, HaloTreeData *&pht) {
    Int_t i,j;
    if (opt.haloidval>0) {
        for (i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
        if (i>=StartSnap && i<EndSnap) {
#endif
            if(pht[i].numhalos>0) for (j=0;j<pht[i].numhalos;j++) pht[i].Halo[j].haloID+=opt.haloidval*(i+opt.snapshotvaloffset)+opt.haloidoffset;
#ifdef USEMPI
        }
#endif
        }
    }
}
//@}

/// \name if particle ids need to be mapped to indices
//@{
void MapPIDStoIndex(Options &opt, HaloTreeData *&pht) {
    cout<<"Mapping PIDS to index "<<endl;
    Int_t i,j,k;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
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

//check to see if ID data compatible with accessing index array allocated
void IDcheck(Options &opt, HaloTreeData *&pht){
    int ierrorflag=0,ierrorsumflag=0;
    for (int i=opt.numsnapshots-1;i>=0;i--) {
#ifdef USEMPI
    if (i>=StartSnap && i<EndSnap) {
#endif
        ierrorflag=0;
        for (Int_t j=0;j<pht[i].numhalos;j++) {
            for (Int_t k=0;k<pht[i].Halo[j].NumberofParticles;k++) 
                if (pht[i].Halo[j].ParticleID[k]<0||pht[i].Halo[j].ParticleID[k]>opt.MaxIDValue) {
                    cout<<"snapshot "<<i<<" particle id out of range "<<pht[i].Halo[j].ParticleID[k]<<" not in [0,"<<opt.MaxIDValue<<")"<<endl;
                    ierrorflag=1;
                    break;
                }
            if (ierrorflag==1) {ierrorsumflag+=1;break;}
        }
#ifdef USEMPI
    }
#endif
    }
#ifdef USEMPI
    int mpierrorflag;
    MPI_Allreduce(&ierrorsumflag,&mpierrorflag,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    ierrorsumflag=mpierrorflag;
#endif
    if (ierrorsumflag>0) {
#ifdef USEMPI
        if (ThisTask==0) 
#endif
        cout<<"Error in particle ids, outside of range. Exiting"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
}

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
