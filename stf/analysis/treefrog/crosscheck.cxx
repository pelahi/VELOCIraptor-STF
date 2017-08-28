
/*! \file crosscheck.cxx
 *  \brief this file contains routines for determining and calculating cross matches
 */

#include "TreeFrog.h"


/// \name cross matching routines
//@{

/*! Determine initial progenitor list based on just merit
  does not guarantee that progenitor list is exclusive.
  The code allows for several types of merit functions and also several
  different ways of using the particle list to calcualte merits.
  By default, all particles belonging to an object are used to find any and all
  matches (if any particle belongs to any object in another catalogue/snapshot)
  This is the default operation (\ref Options.icorematchtype == \ref PARTLISTNOCORE)
  but if particles are arranged in the input catalogue in a meaningful fashion (such
  is the case for VELOCIraptor where particles are arranged in decreasing binding energy)
  then you can use CORE particles set by using fractions of particles and min number
  of particles (\ref Options.particle_frac & \ref Options.min_numpart) combined with
  \ref PARTLISTCORE or \ref PARTLISTCORECORE or \ref PARTLISTCORECOREONLY. If \ref PARTLISTCORE is set, then
  the code uses the core particles to find any and all matches of any of the particles
  belonging to other objects in another catalogue/snapshot. If \re PARTLISTCORECORE
  then matches first done using all particles then for core particle lists in both
  catalogues/snapshots to update the best match. \ref PARTLISTCORECOREONLY only matches cores.
  \todo still need to clean up core matching section and implement \ref PARTLISTCORECOREONLY
*/

ProgenitorData *CrossMatch(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, unsigned int*&pfof2, int &ilistupdated, int istepval, ProgenitorData *refprogen)
{
    long int i,j,k,n,index;
    Int_t numshared;
    Double_t merit;
    int nthreads=1,tid,chunksize;
    long unsigned offset;
    long long hid;
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
private(j,k,tid,pq,numshared,merit,index,offset,np1,np2,pq2,hid)
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
        np1=h1[i].NumberofParticles;
        //if doing core to all matching as set by icorematchtype==PARTLISTCORE
        //then adjust number of particles searched. This does require particles
        //to be meaningfully sorted.
        if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
            np1*=opt.particle_frac;
            if (np1<opt.min_numpart) np1=opt.min_numpart;
            if (np1>h1[i].NumberofParticles) np1=h1[i].NumberofParticles;
        }
        for (j=0;j<np1;j++){
            hid=pfof2[h1[i].ParticleID[j]];
            //correction if use core weighted particles as well.
            if (opt.particle_frac<1 && opt.particle_frac>0 && hid>nhalos2) hid-=nhalos2;
            //store the index in the share and halolist arrays
            index=offset+hid-(long int)1;
            if (hid>0) {
                sharelist[index]+=1;
                //if first time halo has been added, update halolist and increase numshared
                if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
            }
        }
        //now proccess numshared list to remove truly insignificant connections
        for (k=0;k<numshared;k++) {
            j=halolist[offset+k];
            index=offset+j;
            //if sharelist not enough, then move to end and decrease numshared
            np2=h2[j].NumberofParticles;
            //adjust the expecte number of matches if only expecting to match core particles
            if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                np2*=opt.particle_frac;
                if (np2<opt.min_numpart) np2=opt.min_numpart;
                if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
            }
            if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2)) {
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
                np2=h2[j].NumberofParticles;
                //adjust the expecte number of matches if only expecting to match core particles
                if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                    np2*=opt.particle_frac;
                    if (np2<opt.min_numpart) np2=opt.min_numpart;
                    if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
                }
                if (opt.imerittype==NsharedN1N2)
                    merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                else if (opt.imerittype==NsharedN1)
                    merit=(Double_t)sharelist[index]/(Double_t)np1;
                else if (opt.imerittype==Nshared)
                    merit=(Double_t)sharelist[index];
                else if (opt.imerittype==Nsharedcombo)
                    merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
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
            p1[i].istep=istepval;
            delete pq;
        }
        //if the number shared is zero, do nothing
        else {p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;}

        //if matching core to core
        if (opt.icorematchtype==PARTLISTCORECORE && opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
            //calculate number of particles
            np1=(h1[i].NumberofParticles*opt.particle_frac);
            //if halo has enough particle for meaningful most bound particles check, proceed to find merits. Similar to normal merit
            //calculations but have now np2 for calculating merits
            if (np1>=opt.min_numpart) {
                numshared=0;
                for (j=0;j<np1;j++){
                    //only use particles that belong to the core of another group, that is the hid is > nhalos2
                    hid=pfof2[h1[i].ParticleID[j]]-nhalos2;
                    index=offset+hid-(long int)1;
                    if (hid>0) {
                        sharelist[index]+=1;
                        if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
                    }
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
                        if (opt.imerittype==NsharedN1N2)
                            merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        else if (opt.imerittype==NsharedN1)
                            merit=(Double_t)sharelist[index]/(Double_t)np1;
                        else if (opt.imerittype==Nshared)
                            merit=(Double_t)sharelist[index];
                        else if (opt.imerittype==Nsharedcombo)
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
                        p1[i].istep=istepval;
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
        //always search objects that have no progenitors
        if (refprogen[i].NumberofProgenitors==0) {
            num_noprogen+=1;
        }
        // if object has a progenitor, then futher searching depends on
        //what multistep criterion was set
        else {
            //if also applying merit limit then search if does not meet merit threshold
            if (opt.imultsteplinkcrit==MSLCMERIT) {
                if (refprogen[i].Merit[0]<opt.meritlimit) num_noprogen+=1;
            }
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
            //if also applying merit limit then search if does not meet merit threshold
            if (opt.imultsteplinkcrit==MSLCMERIT && refprogen[i].NumberofProgenitors>0) {
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
private(i,j,n,tid,pq,numshared,merit,index,offset,np1,np2,pq2,hid)
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
            np1=h1[i].NumberofParticles;
            if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                np1*=opt.particle_frac;
                if (np1<opt.min_numpart) np1=opt.min_numpart;
                if (np1>h1[i].NumberofParticles) np1=h1[i].NumberofParticles;
            }
            for (j=0;j<np1;j++){
                hid=pfof2[h1[i].ParticleID[j]];
                //correction if use core weighted particles as well.
                if (opt.particle_frac<1 && opt.particle_frac>0 && hid>nhalos2) hid-=nhalos2;
                index=offset+hid-(long int)1;
                if (hid>0) {
                    sharelist[index]+=1;
                    //if first time halo has been added, update halolist and increase numshared
                    if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
                }
            }
            //now proccess numshared list to remove insignificant connections
            for (n=0;n<numshared;n++) {
                j=halolist[offset+n];
                index=offset+j;
                //if sharelist not enough, then move to end and decrease numshared
                np2=h2[j].NumberofParticles;
                //adjust the expecte number of matches if only expecting to match core particles
                if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                    np2*=opt.particle_frac;
                    if (np2<opt.min_numpart) np2=opt.min_numpart;
                    if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
                }
                if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2)) {
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
                    np2=h2[j].NumberofParticles;
                    if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                        np2*=opt.particle_frac;
                        if (np2<opt.min_numpart) np2=opt.min_numpart;
                        if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
                    }
                    if (opt.imerittype==NsharedN1N2)
                        merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                    else if (opt.imerittype==NsharedN1)
                        merit=(Double_t)sharelist[index]/(Double_t)np1;
                    else if (opt.imerittype==Nshared)
                        merit=(Double_t)sharelist[index];
                    else if (opt.imerittype==Nsharedcombo)
                        merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
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
                p1[i].istep=istepval;
                delete pq;
            }
            else {p1[i].ProgenitorList=NULL;p1[i].Merit=NULL;}
            if (opt.icorematchtype==PARTLISTCORECORE && opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
                np1=(h1[i].NumberofParticles*opt.particle_frac);
                if (np1>=opt.min_numpart) {
                    numshared=0;
                    for (j=0;j<np1;j++){
                        //only use particles that belong to the core of another group, that is the hid is > nhalos2
                        hid=pfof2[h1[i].ParticleID[j]]-nhalos2;
                        index=offset+hid-(long int)1;
                        if (hid>0) {
                            sharelist[index]+=1;
                            if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
                        }
                    }
                    for (n=0;n<numshared;n++) {
                        j=halolist[offset+n];
                        index=offset+j;
                        //if sharelist not enough, then move to end and decrease numshared
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
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
                            if (opt.imerittype==NsharedN1N2)
                                merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                            else if (opt.imerittype==NsharedN1)
                                merit=(Double_t)sharelist[index]/(Double_t)np1;
                            else if (opt.imerittype==Nshared)
                                merit=(Double_t)sharelist[index];
                            else if (opt.imerittype==Nsharedcombo)
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
                            p1[i].istep=istepval;
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
    long long hid;
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
private(i,j,k,tid,pq,numshared,merit,index,offset,np1,np2,pq2,hid)
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
        np1=h1[i].NumberofParticles;
        if (opt.icorematchtype>=PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
            np1*=opt.particle_frac;
            if (np1<opt.min_numpart) np1=opt.min_numpart;
            if (np1>h1[i].NumberofParticles) np1=h1[i].NumberofParticles;
        }
        for (j=0;j<np1;j++){
            hid=pfof2[h1[i].ParticleID[j]];
            //correction if use core weighted particles as well.
            if (opt.particle_frac<1 && opt.particle_frac>0 && hid>nhalos2) hid-=nhalos2;
            index=offset+hid-(long int)1;
            if (hid>0) {
                sharelist[index]+=1;
                //if first time halo has been added, update halolist and increase numshared
                if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
            }
        }
        //now calculate merits and sort according to best merit.
        if (numshared>0) {
            //store all viable matches. Merits are calculated here using full number of particles
            np1=h1[i].NumberofParticles;
            pq=new PriorityQueue(numshared);
            for (k=0;k<numshared;k++) {
                j=halolist[offset+k];
                index=offset+j;
                np2=h2[j].NumberofParticles;
                if (opt.imerittype==NsharedN1N2)
                    merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                else if (opt.imerittype==NsharedN1)
                    merit=(Double_t)sharelist[index]/(Double_t)np1;
                else if (opt.imerittype==Nshared)
                    merit=(Double_t)sharelist[index];
                else if (opt.imerittype==Nsharedcombo)
                    merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                pq->Push(j,merit);
                sharelist[index]=0;
            }
            d1[i].NumberofDescendants=numshared;
            d1[i].DescendantList=new long unsigned[numshared];
            d1[i].Merit=new float[numshared];
            d1[i].dtoptype=new unsigned short[numshared];
            //d1[i].nsharedfrac=new float[numshared];
            for (j=0;j<numshared;j++){
                d1[i].DescendantList[j]=pq->TopQueue();
                d1[i].Merit[j]=pq->TopPriority();
                d1[i].dtoptype[j]=1;
                pq->Pop();
            }
            d1[i].istep=istepval;
            delete pq;
        }
        //if the number shared is zero, do nothing
        else {d1[i].DescendantList=NULL;d1[i].Merit=NULL;}

        //if weighted merit function is to be calculated then use the most bound fraction of particles to construct share list
        if (opt.icorematchtype!=PARTLISTNOCORE && opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
            //calculate number of particles
            np1=(h1[i].NumberofParticles*opt.particle_frac);
            if (h1[i].NumberofParticles<opt.min_numpart) np1=h2[j].NumberofParticles;
            else if (np1<opt.min_numpart) np1=opt.min_numpart;
            //if halo has enough particle for meaningful most bound particles check, proceed to find merits. Similar to normal merit
            //calculations but have now np2 for calculating merits
            if (np1>=opt.min_numpart) {
                numshared=0;
                for (j=0;j<np1;j++){
                    hid=pfof2[h1[i].ParticleID[j]]-nhalos2;
                    index=offset+hid-(long int)1;
                    if (hid>0){
                        sharelist[index]+=1;
                        if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
                    }
                }
                if (numshared>0) {
                    //store all viable matches and old ones too
                    pq2=new PriorityQueue(numshared);
                    for (k=0;k<numshared;k++) {
                        j=halolist[offset+k];
                        index=offset+j;
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
                        if (opt.imerittype==NsharedN1N2)
                            merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        else if (opt.imerittype==NsharedN1)
                            merit=(Double_t)sharelist[index]/(Double_t)np1;
                        else if (opt.imerittype==Nshared)
                            merit=(Double_t)sharelist[index];
                        else if (opt.imerittype==Nsharedcombo)
                            merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        pq2->Push(j,merit);
                        sharelist[index]=0;
                    }
                    //update the list, first start with merits
                    for (j=0;j<numshared;j++){
                        for (int k=0;k<d1[i].NumberofDescendants;k++) if (d1[i].DescendantList[k]==pq2->TopQueue())
                        {
                            d1[i].Merit[k]=pq2->TopPriority();
                            d1[i].dtoptype[k]=0;
                            break;
                        }
                        pq2->Pop();
                    }
                    delete pq2;
                    //now redo order based on updated merit
                    pq2=new PriorityQueue(d1[i].NumberofDescendants);
                    pq=new PriorityQueue(d1[i].NumberofDescendants);
                    for (j=0;j<d1[i].NumberofDescendants;j++) {
                        pq->Push(d1[i].dtoptype[j],d1[i].Merit[j]);
                        pq2->Push(d1[i].DescendantList[j],d1[i].Merit[j]);
                    }
                    for (j=0;j<d1[i].NumberofDescendants;j++) {
                        d1[i].DescendantList[j]=pq2->TopQueue();
                        d1[i].Merit[j]=pq2->TopPriority();
                        d1[i].dtoptype[j]=pq->TopQueue();
                        pq2->Pop();
                        pq->Pop();
                    }
                    delete pq;delete pq2;
                    /*
                    for (j=0;j<numshared;j++){
                        d1[i].DescendantList[j]=pq2->TopQueue();
                        d1[i].Merit[j]=pq2->TopPriority();
                        d1[i].dtoptype[j]=0;
                        pq2->Pop();
                    }
                    delete pq2;
                    */
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
        //always search objects that have no progenitors
        if (refdescen[i].NumberofDescendants==0) {
            num_nodescen+=1;
        }
        //if also applying merit limit then search if does not meet merit threshold
        if (opt.imultsteplinkcrit==MSLCPRIMARYPROGEN && refdescen[i].NumberofDescendants>0) {
            if (refdescen[i].dtoptype[0]!=0) num_nodescen+=1;
            else if (refdescen[i].Merit[0]<opt.meritlimit) num_nodescen+=1;
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
            //if also applying merit limit then search if does not meet merit threshold
            if (opt.imultsteplinkcrit==MSLCPRIMARYPROGEN && refdescen[i].NumberofDescendants>0) {
                if (refdescen[i].dtoptype[0]!=0) needdescenlist[num_nodescen++]=i;
                else if (refdescen[i].Merit[0]<opt.meritlimit) needdescenlist[num_nodescen++]=i;
            }
        }

        ntotitems=nhalos2*(long unsigned)nthreads;
        sharelist=new int[ntotitems];
        halolist=new int[ntotitems];
        for (i=0;i<ntotitems;i++)sharelist[i]=0;

#ifdef USEOPENMP
#pragma omp parallel for schedule(dynamic,chunksize) reduction(+:newilistupdated) \
default(shared) \
private(i,j,n,tid,pq,numshared,merit,index,offset,np1,np2,pq2,hid)
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
            np1=h1[i].NumberofParticles;
            if (opt.icorematchtype>=PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                np1*=opt.particle_frac;
                if (np1<opt.min_numpart) np1=opt.min_numpart;
                if (np1>h1[i].NumberofParticles) np1=h1[i].NumberofParticles;
            }
            for (j=0;j<np1;j++){
                hid=pfof2[h1[i].ParticleID[j]];
                //correction if use core weighted particles as well.
                if (opt.particle_frac<1 && opt.particle_frac>0 && hid>nhalos2) hid-=nhalos2;
                index=offset+hid-(long int)1;
                if (hid>0) {
                    sharelist[index]+=1;
                    //if first time halo has been added, update halolist and increase numshared
                    if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
                }
            }
            if (numshared>0) {
                np1=h1[i].NumberofParticles;
                pq=new PriorityQueue(numshared);
                for (n=0;n<numshared;n++){
                    j=halolist[offset+n];
                    index=offset+j;
                    np2=h2[j].NumberofParticles;
                    if (opt.icorematchtype==PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
                        np2*=opt.particle_frac;
                        if (np2<opt.min_numpart) np2=opt.min_numpart;
                        if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
                    }
                    if (opt.imerittype==NsharedN1N2)
                        merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                    else if (opt.imerittype==NsharedN1)
                        merit=(Double_t)sharelist[index]/(Double_t)np1;
                    else if (opt.imerittype==Nshared)
                        merit=(Double_t)sharelist[index];
                    else if (opt.imerittype==Nsharedcombo)
                        merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                    pq->Push(j,merit);
                    sharelist[index]=0;
                }
                //if at this point when looking for progenitors not present in the reference list
                //can check to see if numshared>0 and if so, then reference list will have to be updated
                newilistupdated+=(numshared>0);
                d1[i].NumberofDescendants=numshared;
                d1[i].DescendantList=new long unsigned[numshared];
                d1[i].Merit=new float[numshared];
                d1[i].dtoptype=new unsigned short[numshared];
                //d1[i].nsharedfrac=new float[numshared];
                for (j=0;j<numshared;j++){
                    //d1[i].DescendantList[j]=h2[pq->TopQueue()].haloID;
                    d1[i].DescendantList[j]=pq->TopQueue();
                    d1[i].Merit[j]=pq->TopPriority();
                    d1[i].dtoptype[j]=1;
                    pq->Pop();
                }
                delete pq;
                d1[i].istep=istepval;
            }
            else {d1[i].DescendantList=NULL;d1[i].Merit=NULL;}
            //if weighted merit function is to be calculated then use the most bound fraction of particles to construct share list
            if (opt.icorematchtype!=PARTLISTNOCORE && opt.particle_frac<1 && opt.particle_frac>0 && numshared>0) {
                np1=(h1[i].NumberofParticles*opt.particle_frac);
                if (np1<opt.min_numpart) np1=opt.min_numpart;
                if (np1>h1[i].NumberofParticles) np1=h1[i].NumberofParticles;
                numshared=0;
                for (j=0;j<np1;j++){
                    hid=pfof2[h1[i].ParticleID[j]]-nhalos2;
                    index=offset+hid-(long int)1;
                    if (hid>0) {
                        sharelist[index]+=1;
                        if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
                    }
                }
                if (numshared>0) {
                    //store all viable matches and old ones too
                    pq2=new PriorityQueue(numshared);
                    for (n=0;n<numshared;n++) {
                        j=halolist[offset+n];
                        index=offset+j;
                        np2=h2[j].NumberofParticles*opt.particle_frac;
                        if (h2[j].NumberofParticles<opt.min_numpart) np2=h2[j].NumberofParticles;
                        else if (np2<opt.min_numpart) np2=opt.min_numpart;
                        if (opt.imerittype==NsharedN1N2)
                            merit=(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        else if (opt.imerittype==NsharedN1)
                            merit=(Double_t)sharelist[index]/(Double_t)np1;
                        else if (opt.imerittype==Nshared)
                            merit=(Double_t)sharelist[index];
                        else if (opt.imerittype==Nsharedcombo)
                            merit=(Double_t)sharelist[index]/(Double_t)np1+(Double_t)sharelist[index]*(Double_t)sharelist[index]/(Double_t)np1/(Double_t)np2;
                        pq2->Push(j,merit);
                        sharelist[index]=0;
                    }
                    //update the list, first start with merits
                    for (j=0;j<numshared;j++){
                        for (int k=0;k<d1[i].NumberofDescendants;k++) if (d1[i].DescendantList[k]==pq2->TopQueue())
                        {
                            d1[i].Merit[k]=pq2->TopPriority();
                            d1[i].dtoptype[k]=0;
                            break;
                        }
                        pq2->Pop();
                    }
                    delete pq2;
                    //now redo order based on updated merit
                    pq2=new PriorityQueue(d1[i].NumberofDescendants);
                    pq=new PriorityQueue(d1[i].NumberofDescendants);
                    for (j=0;j<d1[i].NumberofDescendants;j++) {
                        pq->Push(d1[i].dtoptype[j],d1[i].Merit[j]);
                        pq2->Push(d1[i].DescendantList[j],d1[i].Merit[j]);
                    }
                    for (j=0;j<d1[i].NumberofDescendants;j++) {
                        d1[i].DescendantList[j]=pq2->TopQueue();
                        d1[i].Merit[j]=pq2->TopPriority();
                        d1[i].dtoptype[j]=pq->TopQueue();
                        pq2->Pop();
                        pq->Pop();
                    }
                    delete pq;delete pq2;
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
    return d1;
}
//@}

///\name Cleaning cross matches routines
//@{

///ensure cross match is exclusive
void CleanCrossMatch(const int istepval, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, ProgenitorData *&p1)
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
    Double_t *merit=new Double_t[nhalos2];
    //store initial matches
    for (i=0;i<nhalos2;i++) {nh2index[i]=merit[i]=-1;nh2nummatches[i]=0;}//count[i]=0;}
    //ensure that only keep matches at the current time that is being processed and cleaned
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval) {
        for(j=0;j<p1[i].NumberofProgenitors;j++) {
            nh2nummatches[p1[i].ProgenitorList[j]]++;
            if(p1[i].Merit[j]>merit[p1[i].ProgenitorList[j]]){
                merit[p1[i].ProgenitorList[j]]=p1[i].Merit[j];
                nh2index[p1[i].ProgenitorList[j]]=i;
            }
        }
    }

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic,10) nowait
#endif
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval) {
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
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval){
        for(j=0;j<p1[i].NumberofProgenitors;j++) {
            //store nshared
            //p1[i].nsharedfrac[j]=sqrt(p1[i].Merit[j]*((double)h2[p1[i].ProgenitorList[j]].NumberofParticles*(double)h1[i].NumberofParticles))/(double)h1[i].NumberofParticles;
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
///adjust descendant list to store halo ids not halo indices but does not prune the list descendant list.
///\todo really should update to flag/remove all except the truly principle descendant
void CleanCrossMatchDescendant(const int istepval, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, DescendantData *&p1)
{
    Int_t i,j,k;
    int nthreads=1,tid;
    /*
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    int *nh2index=new int[nhalos2];
    int *nh2nummatches=new int[nhalos2];
    Double_t *merit=new Double_t[nhalos2];
    //store initial matches
    for (i=0;i<nhalos2;i++) {nh2index[i]=merit[i]=-1;nh2nummatches[i]=0;}//count[i]=0;}
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval){
        for(j=0;j<p1[i].NumberofDescendants;j++) {
            nh2nummatches[p1[i].DescendantList[j]]++;
            if(p1[i].Merit[j]>merit[p1[i].DescendantList[j]]){
                merit[p1[i].DescendantList[j]]=p1[i].Merit[j];
                nh2index[p1[i].DescendantList[j]]=i;
            }
        }
    }

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic,10) nowait
#endif
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval){
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
    */

    //adjust data to store haloIDS
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic,10) nowait
#endif
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval){
        for(j=0;j<p1[i].NumberofDescendants;j++) {
            //store nshared
            //p1[i].nsharedfrac[j]=sqrt(p1[i].Merit[j]*((double)h2[p1[i].DescendantList[j]].NumberofParticles*(double)h1[i].NumberofParticles))/(double)h1[i].NumberofParticles;
            p1[i].DescendantList[j]=h2[p1[i].DescendantList[j]].haloID;
        }
    }
#ifdef USEOPENMP
}
#endif
    //delete[] nh2index;
    //delete[] nh2nummatches;
    //delete[] merit;
}

//@}

///\name Multistep linking routines for progenitor based searches
//@{
///when linking using multiple snapshots, use to update the progenitor list based on a candidate list built using a single snapshot
///If complex updating done, then will need to update the progenitor based descendant list.
void UpdateRefProgenitors(Options &opt, const Int_t numhalos, ProgenitorData *&pref, ProgenitorData *&ptemp, DescendantDataProgenBased **&pprogendescen, Int_t itime)
{
    if (opt.imultsteplinkcrit==MSLCMISSING) {
        for (Int_t i=0;i<numhalos;i++)
            if (pref[i].NumberofProgenitors==0 && ptemp[i].NumberofProgenitors>0) pref[i]=ptemp[i];
    }
    else if (opt.imultsteplinkcrit==MSLCMERIT) {
        for (Int_t i=0;i<numhalos;i++) {
            //if reference has no links but new links found, copy
            if (pref[i].NumberofProgenitors==0 && ptemp[i].NumberofProgenitors>0) pref[i]=ptemp[i];
            //if reference has no links but new links found also copy as reference did not satisfy merit limit
            else if (pref[i].NumberofProgenitors>0 && ptemp[i].NumberofProgenitors>0) {
                //but only update if new link satisfies merit limit
                if (ptemp[i].Merit[0]>opt.meritlimit) {
                    //first must update the progenitor based descendant list
                    RemoveLinksProgenitorBasedDescendantList(itime, i, pref[i], pprogendescen);
                    //then copy new links
                    pref[i]=ptemp[i];
                }
            }
        }
    }
}

///Construction of descendant list using the candidate progenitor list
void BuildProgenitorBasedDescendantList(Int_t itimeprogen, Int_t itimedescen, Int_t nhalos, ProgenitorData *&pprogen, DescendantDataProgenBased **&pprogendescen, int istep)
{
    //if first pass then store descendants looking at all progenitors
    int did;
    for (Int_t k=0;k<nhalos;k++) {
        for (Int_t iprogen=0;iprogen<pprogen[k].NumberofProgenitors;iprogen++) if (pprogen[k].istep==istep) {
            did=pprogen[k].ProgenitorList[iprogen]-1;//make sure halo descendent index set to start at 0
            pprogendescen[itimedescen][did].haloindex.push_back(k);
            pprogendescen[itimedescen][did].halotemporalindex.push_back(itimeprogen);
            pprogendescen[itimedescen][did].Merit.push_back(pprogen[k].Merit[iprogen]);
            pprogendescen[itimedescen][did].deltat.push_back(istep);
            pprogendescen[itimedescen][did].descentype.push_back(iprogen);
#ifdef USEMPI
            pprogendescen[itimedescen][did].MPITask.push_back(ThisTask);
#endif
            pprogendescen[itimedescen][did].NumberofDescendants++;
        }
    }
}

//removes links of an individual halo
void RemoveLinksProgenitorBasedDescendantList(Int_t itime, Int_t ihaloindex, ProgenitorData &pprogen, DescendantDataProgenBased **&pprogendescen)
{
    //if first pass then store descendants looking at all progenitors
    Int_t itimedescen,did,k=0;
    for (Int_t iprogen=0;iprogen<pprogen.NumberofProgenitors;iprogen++){
        did=pprogen.ProgenitorList[iprogen]-1;//make sure halo descendent index set to start at 0
        itimedescen=itime-pprogen.istep;
        //find where this link exists and then remove it. remove(begin,end, ihaloindex)
        //pprogendescen[itimedescen][did].haloindex.erase(remove(pprogendescen[itimedescen][did].haloindex.begin(),pprogendescen[itimedescen][did].haloindex.end(),ihaloindex),pprogendescen[itimedescen][did].haloindex.end());
        k=0;
        while (k<pprogendescen[itimedescen][did].NumberofDescendants && !(pprogendescen[itimedescen][did].haloindex[k]==ihaloindex && pprogendescen[itimedescen][did].halotemporalindex[k]==itime)) k++;
        pprogendescen[itimedescen][did].haloindex.erase(pprogendescen[itimedescen][did].haloindex.begin()+k);
        pprogendescen[itimedescen][did].halotemporalindex.erase(pprogendescen[itimedescen][did].halotemporalindex.begin()+k);
        pprogendescen[itimedescen][did].Merit.erase(pprogendescen[itimedescen][did].Merit.begin()+k);
        pprogendescen[itimedescen][did].deltat.erase(pprogendescen[itimedescen][did].deltat.begin()+k);
        pprogendescen[itimedescen][did].descentype.erase(pprogendescen[itimedescen][did].descentype.begin()+k);
#ifdef USEMPI
        pprogendescen[itimedescen][did].MPITask.erase(pprogendescen[itimedescen][did].MPITask.begin()+k);
#endif
        pprogendescen[itimedescen][did].NumberofDescendants--;
    }
}

///Cleans the progenitor list to ensure that a given object only has a single descendant using the descendant list built using progenitors. Called if linking across multiple snapshots
void CleanProgenitorsUsingDescendants(Int_t i, HaloTreeData *&pht, DescendantDataProgenBased **&pprogendescen, ProgenitorData **&pprogen, int iopttemporalmerittype)
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
        pprogendescen[i][k].OptimalTemporalMerit(iopttemporalmerittype);
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

///\name Multistep linking routines for descendant based searches
//@{
///similar to \ref UpdateRefProgenitors but for descendants
///\todo need to update how links are removed
void UpdateRefDescendants(Options &opt, const Int_t numhalos, DescendantData *&dref, DescendantData *&dtemp, ProgenitorDataDescenBased **&pdescenprogen, Int_t itime)
{
    if (opt.imultsteplinkcrit==MSLCMISSING) {
        for (Int_t i=0;i<numhalos;i++)
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) dref[i]=dtemp[i];
    }
    else if (opt.imultsteplinkcrit==MSLCPRIMARYPROGEN) {
        for (Int_t i=0;i<numhalos;i++) {
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) dref[i]=dtemp[i];
            else if (dref[i].NumberofDescendants>0 && dtemp[i].NumberofDescendants>0) {
                if (dtemp[i].dtoptype[0]<dref[i].dtoptype[0] || (dref[i].dtoptype[0]>0 && dtemp[i].dtoptype[0]==dref[i].dtoptype[0] && dtemp[i].Merit[0]>dref[i].Merit[0])) {
                    RemoveLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                    //then copy new links
                    dref[i]=dtemp[i];
                    AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                }
            }
        }
    }
}

///similar to construction of descendant list using the candidate progenitor list, but here working in reverse direction
///since objects can have many progenitors but an object should really have only one descendant, the ranking stored in progentype
///is useful for flagging objects which might have complex histories.
void BuildDescendantBasedProgenitorList(Int_t itimedescen, Int_t nhalos, DescendantData *&pdescen, ProgenitorDataDescenBased *&pdescenprogen, int istep)
{
    //if first pass then store progenitors looking at all descendants
    int did;
    for (Int_t k=0;k<nhalos;k++) {
        for (Int_t idescen=0;idescen<pdescen[k].NumberofDescendants;idescen++) if (pdescen[k].istep==istep) {
            did=pdescen[k].DescendantList[idescen]-1;//make sure halo progenitor index set to start at 0
            pdescenprogen[did].haloindex.push_back(k);
            pdescenprogen[did].halotemporalindex.push_back(itimedescen);
            pdescenprogen[did].Merit.push_back(pdescen[k].Merit[idescen]);
            pdescenprogen[did].deltat.push_back(istep);
            pdescenprogen[did].progentype.push_back(idescen);
#ifdef USEMPI
            pdescenprogen[did].MPITask.push_back(ThisTask);
#endif
            pdescenprogen[did].NumberofProgenitors++;
        }
    }
}

///for descedants unless special operation desired, rank the descendant information stored in dtoptype based on the descedant looking backward.
///so if an object had two or more progenitors, the descendant knows its ranking as a progenitor
void UpdateDescendantUsingDescendantBasedProgenitorList(Int_t nhalos,
    DescendantData *&pdescen, ProgenitorDataDescenBased *&pdescenprogen, int istep)
{
    PriorityQueue *pq;
    Int_t nattime;
    Int_t progindex,descenindex,descenprogenindex;
    for (Int_t k=0;k<nhalos;k++) if (pdescenprogen[k].NumberofProgenitors>1)
    {
        //for each halo descandant rank its progenior haloes based on the merit at a given time
        nattime=0;
        for (auto iprogen=0;iprogen<pdescenprogen[k].NumberofProgenitors;iprogen++) nattime+=(pdescenprogen[k].deltat[iprogen]==istep);
        if (nattime<=1) continue;
        //fill the priority queue so that the best merit are first
        pq=new PriorityQueue(nattime);
        for (auto iprogen=0;iprogen<pdescenprogen[k].NumberofProgenitors;iprogen++) if (pdescenprogen[k].deltat[iprogen]==istep)
        {
            pq->Push(iprogen,pdescenprogen[k].Merit[iprogen]);
        }
        nattime=0;
        while(pq->Size()>0)
        {
            progindex=pq->TopQueue();
            descenindex=pdescenprogen[k].haloindex[progindex];
            descenprogenindex=pdescenprogen[k].progentype[progindex];
            pdescen[descenindex].dtoptype[descenprogenindex]=nattime;
            nattime++;
            pq->Pop();
        }
        delete pq;
    }
}

//removes links of an individual halo
void RemoveLinksDescendantBasedProgenitorList(Int_t itime, Int_t ihaloindex, DescendantData &pdescen, ProgenitorDataDescenBased **&pdescenprogen)
{
    //if first pass then store progenitors looking at all descendants
    Int_t itimeprogen,did,k=0;
    for (Int_t idescen=0;idescen<pdescen.NumberofDescendants;idescen++){
        did=pdescen.DescendantList[idescen]-1;//make sure halo descendent index set to start at 0
        itimeprogen=itime+pdescen.istep;
        //find where this link exists and then remove it.
        k=0;
        while (k<pdescenprogen[itimeprogen][did].NumberofProgenitors && !(pdescenprogen[itimeprogen][did].haloindex[k]==ihaloindex && pdescenprogen[itimeprogen][did].halotemporalindex[k]==itime)) k++;
        pdescenprogen[itimeprogen][did].haloindex.erase(pdescenprogen[itimeprogen][did].haloindex.begin()+k);
        pdescenprogen[itimeprogen][did].halotemporalindex.erase(pdescenprogen[itimeprogen][did].halotemporalindex.begin()+k);
        pdescenprogen[itimeprogen][did].Merit.erase(pdescenprogen[itimeprogen][did].Merit.begin()+k);
        pdescenprogen[itimeprogen][did].deltat.erase(pdescenprogen[itimeprogen][did].deltat.begin()+k);
        pdescenprogen[itimeprogen][did].progentype.erase(pdescenprogen[itimeprogen][did].progentype.begin()+k);
#ifdef USEMPI
        pdescenprogen[itimeprogen][did].MPITask.erase(pdescenprogen[itimeprogen][did].MPITask.begin()+k);
#endif
        pdescenprogen[itimeprogen][did].NumberofProgenitors--;
    }
}

//add links of an individual halo
void AddLinksDescendantBasedProgenitorList(Int_t itime, Int_t ihaloindex, DescendantData &pdescen, ProgenitorDataDescenBased **&pdescenprogen)
{
    //if first pass then store progenitors looking at all descendants
    Int_t itimeprogen,did;
    for (Int_t idescen=0;idescen<pdescen.NumberofDescendants;idescen++){
        did=pdescen.DescendantList[idescen]-1;//make sure halo descendent index set to start at 0
        itimeprogen=itime+pdescen.istep;
        pdescenprogen[itimeprogen][did].haloindex.push_back(ihaloindex);
        pdescenprogen[itimeprogen][did].halotemporalindex.push_back(itime);
        pdescenprogen[itimeprogen][did].Merit.push_back(pdescen.Merit[idescen]);
        pdescenprogen[itimeprogen][did].deltat.push_back(pdescen.istep);
        pdescenprogen[itimeprogen][did].progentype.push_back(idescen);
#ifdef USEMPI
        pdescenprogen[itimeprogen][did].MPITask.push_back(ThisTask);
#endif
        pdescenprogen[itimeprogen][did].NumberofProgenitors++;
    }
}

///Update all progenitors of a descandant based on their temporally generalized merit.
void CleanDescendantsUsingProgenitors(Int_t itimeprogen, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdescen, int iopttemporalmerittype)
{
    //if first snapshot can't have any progenitors
    if (itimeprogen==0) return;
    //if no halos before can't have any progenitors
    if (pht[itimeprogen-1].numhalos==0) return;
    PriorityQueue *pq;
    Int_t progindex,descenindex,descentemporalindex,descenprogenindex;
    Double_t generalizedmerit;
    unsigned short rank;
    for (Int_t k=0;k<pht[itimeprogen].numhalos;k++) if (pdescenprogen[itimeprogen][k].NumberofProgenitors>1)
    {
        pq=new PriorityQueue(pdescenprogen[itimeprogen][k].NumberofProgenitors);
        for (auto iprogen=0;iprogen<pdescenprogen[itimeprogen][k].NumberofProgenitors;iprogen++)
        {
            generalizedmerit=pdescenprogen[itimeprogen][k].Merit[iprogen]/pow((Double_t)pdescenprogen[itimeprogen][k].deltat[iprogen],ALPHADELTAT);
            pq->Push(iprogen,generalizedmerit);
        }
        rank=0;
        while(pq->Size()>0)
        {
            progindex=pq->TopQueue();
            descenindex=pdescenprogen[itimeprogen][k].haloindex[progindex];
            descentemporalindex=pdescenprogen[itimeprogen][k].halotemporalindex[progindex];
            descenprogenindex=pdescenprogen[itimeprogen][k].progentype[progindex];
            pdescen[descentemporalindex][descenindex].dtoptype[descenprogenindex]=rank++;
            pq->Pop();
        }
        delete pq;
    }
}

//@}
///Descendant based cleanup routines
//@{
///reranks the data stored in the descendant based on the descendant to progenitor ranking and then the merit.
void RerankDescendants(Options &opt, HaloTreeData *&pht, DescendantData **&pdescen){
#ifndef USEMPI
    int StartSnap=0,EndSnap=opt.numsnapshots;
#endif
    PriorityQueue *pq;
    int mindtop,index;
    Double_t maxmerit;
    unsigned long descen;
    for (int i=0;i<opt.numsnapshots;i++) {
        if (i>=StartSnap && i<EndSnap-1) {
            for (Int_t j=0;j<pht[i].numhalos;j++){
                if (pdescen[i][j].NumberofDescendants>1) {
                    //find the lowest dtop value with highest merit and place it first in the list of descendants
                    mindtop=pdescen[i][j].dtoptype[0];
                    maxmerit=pdescen[i][j].Merit[0];
                    descen=pdescen[i][j].DescendantList[0];
                    index=0;
                    //if a zero rank is first, no need to do anything
                    if (mindtop==0) continue;
                    for (int k=1;k<pdescen[i][j].NumberofDescendants;k++) {
                        if (pdescen[i][j].dtoptype[k]<mindtop || (pdescen[i][j].dtoptype[k]==mindtop&&pdescen[i][j].Merit[k]>maxmerit)) {
                            mindtop=pdescen[i][j].dtoptype[k];
                            maxmerit=pdescen[i][j].Merit[k];
                            descen=pdescen[i][j].DescendantList[k];
                            index=k;
                        }
                    }
                    //swap data
                    if (index!=0) {
                        pdescen[i][j].DescendantList[index]=pdescen[i][j].DescendantList[0];
                        pdescen[i][j].Merit[index]=pdescen[i][j].Merit[0];
                        pdescen[i][j].dtoptype[index]=pdescen[i][j].dtoptype[0];
                        pdescen[i][j].DescendantList[0]=descen;
                        pdescen[i][j].Merit[0]=maxmerit;
                        pdescen[i][j].dtoptype[0]=mindtop;
                    }
                }
            }
        }
    }

}
//@}
//@}
