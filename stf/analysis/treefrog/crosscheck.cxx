
/*! \file crosscheck.cxx
 *  \brief this file contains routines for determining and calculating cross matches
 */

#include "TreeFrog.h"


/// \name cross matching routines
//@{

///Calculate the merit between two haloes that have been matched
Double_t CalculateMerit(Options &opt, UInt_t n1, UInt_t n2, UInt_t nsh, HaloData &h1, HaloData &h2, UInt_t hindex=0, unsigned int *sharepartlist=NULL, Double_t *rankingsum=NULL)
{
    Double_t merit;
    Double_t ranksum,norm;
    if (opt.imerittype==MERITNsharedN1N2)
        merit=(Double_t)nsh*(Double_t)nsh/(Double_t)n1/(Double_t)n2;
    else if (opt.imerittype==MERITNsharedN1)
        merit=(Double_t)nsh/(Double_t)n1;
    else if (opt.imerittype==MERITNshared)
        merit=(Double_t)nsh;
    else if (opt.imerittype==MERITNsharedcombo)
        merit=(Double_t)nsh/(Double_t)n1*(1.0+(Double_t)nsh/(Double_t)n2);
    else if (opt.imerittype==MERITRankWeighted) {
        //this ranking is based on Poole+ 2017 where we rank most bound particles more
        //assumes input particle list is meaningful (either boundness ranked or radially for example)
        ranksum=0;
        for (auto i=0;i<n1;i++) if (sharepartlist[i]==hindex) ranksum+=1.0/(i+1.0);
        //normalize to the optimal value for nsh=n2, all particles in the descendant
        norm=0.5772156649+log((Double_t)n2);
        merit=ranksum/norm;
        //and multiply this to standard Nsh^2/N1/N2 merit to correct for small objects
        //being lost in bigger object and being deemed main progenitor
        merit*=(Double_t)nsh*(Double_t)nsh/(Double_t)n1/(Double_t)n2;
        merit=sqrt(merit);
    }
    else if (opt.imerittype==MERITRankWeightedBoth) {
        //like above but ranking both ways, that is merit is combination of rankings in a and b
        ranksum=0;
        for (auto i=0;i<n1;i++) if (sharepartlist[i]==hindex) ranksum+=1.0/(i+1.0);
        //normalize to the optimal value for nsh=n2, all particles in the descendant
        norm=0.5772156649+log((Double_t)n1);
        merit=ranksum/norm;
        ranksum=rankingsum[hindex-1];
        //reset value after use
        rankingsum[hindex-1]=0;
        norm=0.5772156649+log((Double_t)n2);
        merit*=ranksum/norm;
        merit*=(Double_t)nsh*(Double_t)nsh/(Double_t)n1/(Double_t)n2;
        merit=sqrt(merit);
    }
    return merit;
}

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
    long int i,j,k,n;
    int nthreads=1,tid,maxnthreads=1,chunksize;
    long unsigned offset;
    //temp variable to store the openmp reduction value for ilistupdated
    int newilistupdated;
    //setup openmp environment
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
    }
    //set the active number of threads based on breaking up computation in large enough chunks
    nthreads=min((int)ceil((double)nhalos1/(double)OMPCHUNKSIZE),maxnthreads);
    omp_set_num_threads(nthreads);
    chunksize=OMPCHUNKSIZE;
#endif

    ProgenitorData *p1=new ProgenitorData[nhalos1];
    unsigned int *sharelist,*halolist;
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
    sharelist=new unsigned int[ntotitems];
    halolist=new unsigned int[ntotitems];
    //to store haloes that share links and the halo index of those shared haloes
    for (i=0;i<ntotitems;i++)sharelist[i]=0;

    offset=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid,offset)
{
        //initialize variables
        tid=omp_get_thread_num();
        offset=((long unsigned)tid)*nhalos2;
#pragma omp for schedule(dynamic,chunksize) nowait
#endif
    for (i=0;i<nhalos1;i++){
        CrossMatchProgenitorIndividual(opt, i, nhalos1, nhalos2, h1, h2,  pfof2, istepval, p1, sharelist, halolist, offset);
    }
#ifdef USEOPENMP
}
#endif
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
        //set OpenMP environment
#ifdef USEOPENMP
        //set the active number of threads based on breaking up computation in large enough chunks
        nthreads=min((int)ceil((double)num_noprogen/(double)OMPCHUNKSIZE),maxnthreads);
        omp_set_num_threads(nthreads);
#endif
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
        sharelist=new unsigned int[ntotitems];
        halolist=new unsigned int[ntotitems];
        for (i=0;i<ntotitems;i++)sharelist[i]=0;

        offset=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,k,tid,offset)
{
        //initialize variables
        tid=omp_get_thread_num();
        offset=((long unsigned)tid)*nhalos2;
#pragma omp parallel for schedule(dynamic,chunksize) reduction(+:newilistupdated)
#endif
        for (k=0;k<num_noprogen;k++){
            i=needprogenlist[k];
            newilistupdated+=CrossMatchProgenitorIndividual(opt, i, nhalos1, nhalos2, h1, h2,  pfof2, istepval, p1, sharelist, halolist, offset);
        }
#ifdef USEOPENMP
}
#endif
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

    //reset number of threads back to maximum
#ifdef USEOPENMP
#pragma omp parallel
        omp_set_num_threads(maxnthreads);
#endif
    return p1;
}

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
    )
{
    long int j,k,index;
    Int_t numshared;
    Double_t merit;
    long long hid;
    PriorityQueue *pq,*pq2;
    UInt_t np1,np2;
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
        if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2) && sharelist[index]<opt.mlsig*sqrt((Double_t)np1)) {
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
            merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j]);
            pq->Push(j,merit);
            sharelist[index]=0;
        }
        p1[i].NumberofProgenitors=numshared;
        p1[i].ProgenitorList=new long unsigned[numshared];
        p1[i].Merit=new float[numshared];
        p1[i].nsharedfrac=new float[numshared];
        for (j=0;j<numshared;j++){
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
                if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2) && sharelist[index]<opt.mlsig*sqrt((Double_t)np1)){
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
                    merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j]);
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
    return (numshared>0);
}

///effectively the same code as \ref CrossMatch but allows for the possibility of matching descendant/child nodes using a different merit function
///here the lists are different as input order is reverse of that used in \ref CrossMatch
DescendantData *CrossMatchDescendant(Options &opt, const long unsigned nhalos1, const long unsigned nhalos2,
    HaloData *&h1, HaloData *&h2, unsigned int *&pfof2, int &ilistupdated, int istepval, unsigned int *pranking2,
    DescendantData *refdescen)
{
    long int i,j,k,n;
    int nthreads=1,tid,maxnthreads=1,chunksize;
    long unsigned offset, offset2;
    //temp variable to store the openmp reduction value for ilistupdated
    int newilistupdated;
    int initdtopval;
    //setup openmp environment
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    if (omp_get_thread_num()==0) maxnthreads=omp_get_num_threads();
    }
    //set the active number of threads based on breaking up computation in large enough chunks
    nthreads=min((int)ceil((double)nhalos1/(double)OMPCHUNKSIZE),maxnthreads);
    omp_set_num_threads(nthreads);
    chunksize=OMPCHUNKSIZE;
#endif

    DescendantData *d1=new DescendantData[nhalos1];
    unsigned int *sharelist, *halolist, *sharepartlist=NULL;
    Double_t *rankingsum=NULL;
    UInt_t nbiggest;
    long unsigned num_nodescen, ntotitems;
    long unsigned *needdescenlist;
    if (refdescen==NULL) ilistupdated=1;
    else ilistupdated=0;
    newilistupdated=ilistupdated;
    ntotitems=nhalos2*(long unsigned)nthreads;

    //if a core-core matching follows core-all matching then mark core-all matches
    //as worse initial rank of 1, otherwise use zero
    initdtopval=(opt.icorematchtype!=PARTLISTNOCORE && opt.particle_frac<1 && opt.particle_frac>0);

    if (nhalos2>0){

    //if need to store ranking of shared particles, just identify largest halo and allocate mem to store that object
    if (opt.imerittype==MERITRankWeighted||opt.imerittype==MERITRankWeightedBoth) {
        nbiggest=h1[0].NumberofParticles;
        for (i=1;i<nhalos1;i++) if (nbiggest<h1[i].NumberofParticles) nbiggest=h1[i].NumberofParticles;
        sharepartlist=new unsigned int[nbiggest*(long unsigned)nthreads];
    }
    if (opt.imerittype==MERITRankWeightedBoth) {
        rankingsum=new Double_t[ntotitems];
        for (j=0;j<ntotitems;j++) rankingsum[j]=0;
    }

    if (refdescen==NULL) {
    sharelist=new unsigned int[ntotitems];
    halolist=new unsigned int[ntotitems];
    //to store haloes that share links and the halo index of those shared haloes
    for (i=0;i<ntotitems;i++)sharelist[i]=0;
    offset=offset2=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid,offset,offset2)
{
        //initialize variables
        tid=omp_get_thread_num();
        offset=((long unsigned)tid)*nhalos2;
        offset2=((long unsigned)tid)*nbiggest;
#pragma omp for schedule(dynamic,chunksize) nowait
#endif
    for (i=0;i<nhalos1;i++){
        CrossMatchDescendantIndividual(opt, i, nhalos1, nhalos2, h1, h2, pfof2, istepval, initdtopval, d1,
            sharelist, halolist, offset, offset2, sharepartlist, pranking2, rankingsum);
    }
#ifdef USEOPENMP
}
#endif
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
        //if also applying merit limit then search if does not meet merit threshold or not ideal progenitor
        if (opt.imultsteplinkcrit==MSLCPRIMARYPROGEN && refdescen[i].NumberofDescendants>0) {
            if (refdescen[i].dtoptype[0]!=0) num_nodescen+=1;
        }
        else if (opt.imultsteplinkcrit==MSLCMERITPRIMARYPROGEN && refdescen[i].NumberofDescendants>0) {
            if (refdescen[i].dtoptype[0]!=0) num_nodescen+=1;
            else if (refdescen[i].Merit[0]<opt.meritlimit) num_nodescen+=1;
        }
    }
    //only allocate memory and process list if there are any haloes needing to be searched
    if (num_nodescen>0) {
        //set OpenMP environment
#ifdef USEOPENMP
        //set the active number of threads based on breaking up computation in large enough chunks
        nthreads=min((int)ceil((double)num_nodescen/(double)OMPCHUNKSIZE),maxnthreads);
        omp_set_num_threads(nthreads);
#endif
        needdescenlist=new long unsigned[num_nodescen];
        num_nodescen=0;
        for (i=0;i<nhalos1;i++){
            if (refdescen[i].NumberofDescendants==0){
                needdescenlist[num_nodescen++]=i;
            }
            if (opt.imultsteplinkcrit==MSLCPRIMARYPROGEN && refdescen[i].NumberofDescendants>0) {
                if (refdescen[i].dtoptype[0]!=0) needdescenlist[num_nodescen++]=i;
            }
            else if (opt.imultsteplinkcrit==MSLCMERITPRIMARYPROGEN && refdescen[i].NumberofDescendants>0) {
                if (refdescen[i].dtoptype[0]!=0) needdescenlist[num_nodescen++]=i;
                else if (refdescen[i].Merit[0]<opt.meritlimit) needdescenlist[num_nodescen++]=i;
            }
        }

        sharelist=new unsigned int[ntotitems];
        halolist=new unsigned int[ntotitems];
        for (i=0;i<ntotitems;i++)sharelist[i]=0;

        offset=offset2=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,k,tid,offset,offset2)
{
        //initialize variables
        tid=omp_get_thread_num();
        offset=((long unsigned)tid)*nhalos2;
        offset2=((long unsigned)tid)*nbiggest;
#pragma omp parallel for schedule(dynamic,chunksize) reduction(+:newilistupdated)
#endif
        for (k=0;k<num_nodescen;k++){
            i=needdescenlist[k];
            newilistupdated+=CrossMatchDescendantIndividual(opt, i, nhalos1, nhalos2, h1, h2, pfof2, istepval, initdtopval, d1,
                sharelist, halolist, offset, offset2, sharepartlist, pranking2, rankingsum);
        }
#ifdef USEOPENMP
}
#endif
        delete[] sharelist;
        delete[] halolist;
        delete[] needdescenlist;
    }
    //end of if statement for progenitors more than one snap ago if any current haloes need to be searched
    ilistupdated=newilistupdated;
    }
    //end of progenitors more than one snap ago

    //free some memory
    if (opt.imerittype==MERITRankWeighted||opt.imerittype==MERITRankWeightedBoth) delete[] sharepartlist;
    if (opt.imerittype==MERITRankWeightedBoth) delete[] rankingsum;
    }
    //end of more than zero halos to produce links to
    else {
        for (i=0;i<nhalos1;i++){
            d1[i].NumberofDescendants=0;
            d1[i].DescendantList=NULL;d1[i].Merit=NULL;
        }
    }
    //reset number of threads back to maximum
#ifdef USEOPENMP
#pragma omp parallel
        omp_set_num_threads(maxnthreads);
#endif
    return d1;
}

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
    )
{
    long int j,k,index;
    Int_t numshared;
    Double_t merit;
    long long hid;
    PriorityQueue *pq,*pq2;
    UInt_t np1,np2;
    numshared=0;
    //go through halos particle list to see if these particles belong to another halo
    //at a different time/in a different catalog
    np1=h1[i].NumberofParticles;
    if (opt.icorematchtype>=PARTLISTCORE && opt.particle_frac<1 && opt.particle_frac>0) {
        np1*=opt.particle_frac;
        if (np1<opt.min_numpart) np1=opt.min_numpart;
        if (np1>h1[i].NumberofParticles) np1=h1[i].NumberofParticles;
    }
    //initialize the shared particle list if desired
    if (opt.imerittype==MERITRankWeighted||opt.imerittype==MERITRankWeightedBoth) for (j=0;j<np1;j++) sharepartlist[j+offset2]=0;
    //loop over particles
    for (j=0;j<np1;j++){
        hid=pfof2[h1[i].ParticleID[j]];
        //correction if use core weighted particles as well.
        if (opt.icorematchtype!=PARTLISTNOCORE &&  opt.particle_frac<1 && opt.particle_frac>0 && hid>nhalos2) hid-=nhalos2;
        index=offset+hid-(long int)1;
        if (hid>0) {
            sharelist[index]+=1;
            //if first time halo has been added, update halolist and increase numshared
            if (sharelist[index]==1) halolist[offset+numshared++]=hid-1;
            //if storing the ranking of shared particles
            if (opt.imerittype==MERITRankWeighted||opt.imerittype==MERITRankWeightedBoth) sharepartlist[j+offset2]=hid;
            //if need ranking the other way, calculate that as well
            if (opt.imerittype==MERITRankWeightedBoth) rankingsum[hid-1+offset]+=1.0/(1.0+pranking2[h1[i].ParticleID[j]]);
        }
    }
    //now proccess numshared list to remove insignificant connections
    for (k=0;k<numshared;k++) {
        j=halolist[offset+k];
        index=offset+j;
        //if sharelist not enough, then move to end and decrease numshared
        np2=h2[j].NumberofParticles;
        if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2) && sharelist[index]<opt.mlsig*sqrt((Double_t)np1)) {
            if (opt.imerittype==MERITRankWeightedBoth) rankingsum[j+offset]=0;
            halolist[offset+k]=halolist[offset+numshared-1];
            sharelist[index]=0;
            k--;
            numshared--;
        }
    }
    //now calculate merits and sort according to best merit.
    if (numshared>0) {
        ////store all viable matches. Merits are calculated here using full number of particles
        //np1=h1[i].NumberofParticles;
        pq=new PriorityQueue(numshared);
        for (k=0;k<numshared;k++) {
            j=halolist[offset+k];
            index=offset+j;
            np2=h2[j].NumberofParticles;
            //for openmp compatability, have several if statements
            if (opt.imerittype==MERITRankWeighted) merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j],j+1,&sharepartlist[offset2]);
            else if (opt.imerittype==MERITRankWeightedBoth) merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j],j+1,&sharepartlist[offset2],&rankingsum[offset]);
            else merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j],j+1);
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
            d1[i].dtoptype[j]=initdtopval;
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
        if (np1<opt.min_numpart) np1=opt.min_numpart;
        if (h1[i].NumberofParticles<np1) np1=h1[i].NumberofParticles;
        if (opt.imerittype==MERITRankWeighted||opt.imerittype==MERITRankWeightedBoth) for (j=0;j<np1;j++) sharepartlist[j+offset2]=0;
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
                    if (opt.imerittype==MERITRankWeighted||opt.imerittype==MERITRankWeightedBoth) sharepartlist[j+offset2]=hid;
                    //if need ranking the other way, calculate that as well
                    if (opt.imerittype==MERITRankWeightedBoth) rankingsum[hid-1+offset]+=1.0/(1.0+pranking2[h1[i].ParticleID[j]]);
                }
            }
            //now proccess numshared list to remove insignificant connections
            for (k=0;k<numshared;k++) {
                j=halolist[offset+k];
                index=offset+j;
                //if sharelist not enough, then move to end and decrease numshared
                np2=h2[j].NumberofParticles*opt.particle_frac;
                if (np2<opt.min_numpart) np2=opt.min_numpart;
                if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
                if(sharelist[index]<opt.mlsig*sqrt((Double_t)np2) && sharelist[index]<opt.mlsig*sqrt((Double_t)np1)) {
                    if (opt.imerittype==MERITRankWeightedBoth) rankingsum[j+offset]=0;
                    halolist[offset+k]=halolist[offset+numshared-1];
                    sharelist[index]=0;
                    k--;
                    numshared--;
                }
            }
            if (numshared>0) {
                //store all viable matches and old ones too
                pq2=new PriorityQueue(numshared);
                for (k=0;k<numshared;k++) {
                    j=halolist[offset+k];
                    index=offset+j;
                    np2=h2[j].NumberofParticles*opt.particle_frac;
                    if (np2<opt.min_numpart) np2=opt.min_numpart;
                    if (np2>h2[j].NumberofParticles) np2=h2[j].NumberofParticles;
                    if (opt.imerittype==MERITRankWeighted) merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j],j+1,&sharepartlist[offset2]);
                    else if (opt.imerittype==MERITRankWeightedBoth) merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j],j+1,&sharepartlist[offset2],&rankingsum[offset]);
                    else merit=CalculateMerit(opt,np1,np2,sharelist[index],h1[i],h2[j],j+1);
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
            //end of numshared>0
        }
        //end of halo has enough particles for weighted (bound fraction) merit adjustment
    }
    //end of weighted (based on most bound particles) merit
    return (numshared>0);
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
void UpdateDescendantIndexing(const int istepval, const long unsigned nhalos1, const long unsigned nhalos2, HaloData *&h1, HaloData *&h2, DescendantData *&p1)
{
    Int_t i,j,k;
    int nthreads=1,tid;
    //adjust data to store haloIDS
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k)
{
#pragma omp for schedule(dynamic,10) nowait
#endif
    for (i=0;i<nhalos1;i++) if (p1[i].istep==istepval){
        for(j=0;j<p1[i].NumberofDescendants;j++) {
            p1[i].DescendantList[j]=h2[p1[i].DescendantList[j]].haloID;
        }
    }
#ifdef USEOPENMP
}
#endif
}

///now adjust the rankings of the progenitor descendant matches. Here the logic is applied so that if
///objects are considered the best match for multiple descendants, the rank is adjusted and the
///ranks of the progenitors of these secondary descendants are set to 0
///also then searches for any objects that have more than one possible descendant but none are rank 0.
///these are checked to see if a zero rank descendant can be identified.
void CleanCrossMatchDescendant(Options &opt, Int_t itime, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdescen)
{
    Int_t i,j,k;
    int nthreads=1,tid;
    int iflag;
    int itimedescen, itimeprogen;
    unsigned long did;
    Int_t progindex,descenindex,descentemporalindex,descenprogenindex;
    PriorityQueue *pq;
    unsigned short rank;
    Double_t generalizedmerit, meritratio;
    int numcurprogen;
    int numcorrected=0;
    vector<int> descenindexvec;

    //the idea here is to adjust rankings of descendant to progenitor or remove a connection completely
    //if an object has two our more descendants of rank 0
    for (Int_t k=0;k<pht[itime].numhalos;k++) if (pdescen[itime][k].NumberofDescendants>1)
    {
        //if object has very low merit, ignore
        if (pdescen[itime][k].Merit[0]<opt.meritlimit) continue;
        //find first number of rank 0 descendants
        for (auto idescen=0;idescen<pdescen[itime][k].NumberofDescendants;idescen++)
        {
            if (pdescen[itime][k].dtoptype[idescen]==0 && pdescen[itime][k].Merit[idescen]>opt.meritlimit) {
                descenindexvec.push_back(idescen);
            }
        }
        if (descenindexvec.size()<=1) {
            descenindexvec.clear();
            continue;
        }
        //now search all other objects to see if they have descendants that link to other progenitors
        //if an object has multiple descendants of similar rank, may want to flag it as poor, reset ranks

        for (auto idescen=1;idescen<descenindexvec.size();idescen++) {
            itimedescen=itime+pdescen[itime][k].istep;
            did=pdescen[itime][k].DescendantList[descenindexvec[idescen]]-1;

            //if object has other possible matches, adjust the rank of this descendant
            pdescen[itime][k].dtoptype[descenindexvec[idescen]]=1;
            //find where this information is stored in its descendant and adjust it as well.
            for (auto iprogen=0;iprogen<pdescenprogen[itimedescen][did].NumberofProgenitors;iprogen++) {
                descenindex=pdescenprogen[itimedescen][did].haloindex[iprogen];
                descentemporalindex=pdescenprogen[itimedescen][did].halotemporalindex[iprogen];
                descenprogenindex=pdescenprogen[itimedescen][did].progenindex[iprogen];
                //found object in ProgenitorDataDescenBased that must be changed
                if (descenindex==k && descentemporalindex==itime) {
                    pdescenprogen[itimedescen][did].dtoptype[iprogen]=1;
                    break;
                }
            }

            //if descendant itself does not have any other progenitors do nothing
            if (pdescenprogen[itimedescen][did].NumberofProgenitors<=1) continue;
            //if object does not have another useful progenitor do nothing other than note this object
            iflag=0;
            for (auto iprogen=0;iprogen<pdescenprogen[itimedescen][did].NumberofProgenitors;iprogen++) iflag+=(pdescenprogen[itimedescen][did].Merit[iprogen]>opt.meritlimit);
            if (iflag<=1) continue;


            //otherwise see if the ranking can be swapped.
            //check other possible progenitors of this descendant, compare the merits
            iflag=0;
            for (auto iprogen=0;iprogen<pdescenprogen[itimedescen][did].NumberofProgenitors;iprogen++) {
                descenindex=pdescenprogen[itimedescen][did].haloindex[iprogen];
                descentemporalindex=pdescenprogen[itimedescen][did].halotemporalindex[iprogen];
                descenprogenindex=pdescenprogen[itimedescen][did].progenindex[iprogen];
                //found the object whose rank can be changed. if merit is reasonable, change ranking
                if (!(descenindex==k && descentemporalindex==itime) && pdescenprogen[itimedescen][did].dtoptype[iprogen]==1)
                {
                    meritratio=pdescenprogen[itimedescen][did].Merit[iprogen]/pdescen[itime][k].Merit[descenindexvec[idescen]];
                    if (meritratio<opt.meritratiolimit && 1.0/meritratio<opt.meritratiolimit) {
                        iflag=1;
                        pdescenprogen[itimedescen][did].dtoptype[iprogen]=0;
#ifdef USEMPI
                        if (pdescenprogen[itimedescen][did].MPITask[iprogen]==ThisTask)
#endif
                        pdescen[descentemporalindex][descenindex].dtoptype[descenprogenindex]=0;
                    }
                    break;
                }
            }
            if (iflag==0) continue;
            numcorrected++;
        }
        descenindexvec.clear();
    }
    if (opt.iverbose>=2) cout<<"Number of corrected haloes "<<numcorrected<<endl;
}

void CleanDescendantsForMissingProgenitors(Options &opt, Int_t itime, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdescen)
{
    Int_t i,j,k;
    int nthreads=1,tid;
    int iflag,irank,index,index2,indexprime;
    int itimedescen, itimeprogen,itimedescenprog;
    Int_t descenindex,descentemporalindex,descenprogenindex,descenprogenprimeindex;
    Int_t descenindex2,descentemporalindex2,descenprogenindex2;
    Int_t descenprogindex,descenprogtemporalindex;
    Double_t merit1,merit2,meritprime;
    int numcorrected=0, numall=0;
    //and also check to see if there are any objects with no zero rank progenitors.
    //if this object shares a progenitor with another object that has another reasonable rank progenitor
    //change the ranking
    for (Int_t k=0;k<pht[itime].numhalos;k++) if (pdescenprogen[itime][k].NumberofProgenitors>0)
    {
        //see if object has no zero rank progenitor find best ranked progenitor of merit ordered values
        index=0;
        iflag=(pdescenprogen[itime][k].dtoptype[index]==0);
        irank=pdescenprogen[itime][k].dtoptype[index];
        for (auto iprogen=1;iprogen<pdescenprogen[itime][k].NumberofProgenitors;iprogen++)
        {
            if (irank<pdescenprogen[itime][k].dtoptype[iprogen]) {
                irank=pdescenprogen[itime][k].dtoptype[iprogen];
                index=iprogen;
            }
            if (pdescenprogen[itime][k].dtoptype[iprogen]==0) iflag++;
        }
        if (iflag) continue;
        numall++;
        //now have halo of interest, see if its lowest (best) rank progenitor has more than one descendant
        merit1=pdescenprogen[itime][k].Merit[index];
        //check to see if merit is viable enough to warrant adjusting rankings
        if (merit1<opt.meritlimit) continue;

        //get the progenitor and see if it has other descendants
        descenindex=pdescenprogen[itime][k].haloindex[index];
        descentemporalindex=pdescenprogen[itime][k].halotemporalindex[index];
        descenprogenindex=pdescenprogen[itime][k].progenindex[index];
        if (pdescen[descentemporalindex][descenindex].NumberofDescendants<=1) continue;

        //search the other progenitor's descendants for a zero rank progenitor
        iflag=1;
        for (auto idescen=0;idescen<pdescen[descentemporalindex][descenindex].NumberofDescendants;idescen++) {
            if (pdescen[descentemporalindex][descenindex].dtoptype[idescen]==0) {
                iflag=0;
                descenprogenprimeindex=idescen;
                break;
            }
        }
        if (iflag) continue;

        //look at other descendant to see if it has more than one progenitor
        descenprogindex=pdescen[descentemporalindex][descenindex].DescendantList[descenprogenprimeindex]-1;
        itimedescenprog=pdescen[descentemporalindex][descenindex].istep+descentemporalindex;
        if (pdescenprogen[itimedescenprog][descenprogindex].NumberofProgenitors<=1) continue;

        //now have best rank, lets check merits
        meritprime=pdescen[descentemporalindex][descenindex].Merit[descenprogenprimeindex];
        //if merit is not within some factor then stop
        if (merit1/meritprime>opt.meritratiolimit || meritprime/merit1>opt.meritratiolimit) continue;

        //find where progenitor is in the progenitor  list
        for (auto iprogen=0;iprogen<pdescenprogen[itimedescenprog][descenprogindex].NumberofProgenitors;iprogen++) {
            if (pdescenprogen[itimedescenprog][descenprogindex].haloindex[iprogen]==descenindex && pdescenprogen[itimedescenprog][descenprogindex].halotemporalindex[iprogen]==descentemporalindex) {
                indexprime=iprogen;
                break;
            }
        }
        for (auto iprogen=0;iprogen<pdescenprogen[itimedescenprog][descenprogindex].NumberofProgenitors;iprogen++) {
                int descenindex3=pdescenprogen[itimedescenprog][descenprogindex].haloindex[iprogen];
                int descentemporalindex3=pdescenprogen[itimedescenprog][descenprogindex].halotemporalindex[iprogen];
                int descenprogenindex3=pdescenprogen[itimedescenprog][descenprogindex].progenindex[iprogen];
        }

        //find the index of the other possible progenitor
        for (auto iprogen=0;iprogen<pdescenprogen[itimedescenprog][descenprogindex].NumberofProgenitors;iprogen++) {
            if (pdescenprogen[itimedescenprog][descenprogindex].dtoptype[iprogen]==1) {
                descenindex2=pdescenprogen[itimedescenprog][descenprogindex].haloindex[iprogen];
                descentemporalindex2=pdescenprogen[itimedescenprog][descenprogindex].halotemporalindex[iprogen];
                descenprogenindex2=pdescenprogen[itimedescenprog][descenprogindex].progenindex[iprogen];
                index2=iprogen;
                break;
            }
        }

        merit2=pdescen[descentemporalindex2][descenindex2].Merit[descenprogenindex2];
        //also check other merit to see if swap is allowed
        if (merit2<opt.meritlimit) continue;
        if (merit2/meritprime>opt.meritratiolimit || meritprime/merit2>opt.meritratiolimit) continue;

        //first the rank 1 progenitor of the secondary descendant now becomes the main progenitor
        pdescen[descentemporalindex2][descenindex2].dtoptype[descenprogenindex2]=0;
        pdescenprogen[itimedescenprog][descenprogindex].dtoptype[index2]=0;
        //then we adjust the rank of the previous main progenitor
        pdescen[descentemporalindex][descenindex].dtoptype[descenprogenprimeindex]=1;
        pdescenprogen[itimedescenprog][descenprogindex].dtoptype[indexprime]=1;
        //and finally adjust the original progenitor
        pdescen[descentemporalindex][descenindex].dtoptype[descenprogenindex]=0;
        pdescenprogen[itime][k].dtoptype[index]=0;
        numcorrected++;
    }
    if (opt.iverbose>=2) cout<<"Number of corrected haloes that did not have primary progenitors "<<numcorrected<<" with all possibilities being "<<numall<<endl;
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
#ifdef USEMPI
        //for mpi communication need to store if anything is removed
        pprogendescen[itimedescen][did].removalhaloindex.push_back(pprogendescen[itimedescen][did].haloindex[k]);
        pprogendescen[itimedescen][did].removalhalotemporalindex.push_back(pprogendescen[itimedescen][did].halotemporalindex[k]);
#endif
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
void UpdateRefDescendants(Options &opt, const Int_t numhalos, DescendantData *&dref, DescendantData *&dtemp, ProgenitorDataDescenBased **&pdescenprogen, Int_t itime)
{
    int iprogen,idescen;
    int itimeprogen,itimedescen;
    int iflag;
    //remove all connections present in dtemp after ranking, and only add them if required
    for (Int_t i=0;i<numhalos;i++) if (dtemp[i].NumberofDescendants>0) {
        RemoveLinksDescendantBasedProgenitorList(itime, i, dtemp[i], pdescenprogen);
    }
    //only add new descendant if reference has no descendant
    if (opt.imultsteplinkcrit==MSLCMISSING) {
        for (Int_t i=0;i<numhalos;i++) {
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) {
                dref[i]=dtemp[i];
                AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
            }
        }
    }
    //also add if newly identified descendant link is 0 (primary) and previous was not.
    //note that due to links ony being ranked at a given time, also need to check
    //that descendant does not have a 0 rank progenitor already.
    else if (opt.imultsteplinkcrit==MSLCPRIMARYPROGEN) {
        for (Int_t i=0;i<numhalos;i++) {
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) {
                dref[i]=dtemp[i];
                AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
            }
            else if (dref[i].NumberofDescendants>0 && dtemp[i].NumberofDescendants>0) {
                if (dtemp[i].dtoptype[0]==0 && dref[i].dtoptype[0]!=0) {
                    //examine identified progenitor and see if it already has a 0 rank link
                    idescen=dtemp[i].DescendantList[0]-1;
                    itimedescen=dtemp[i].istep+itime;
                    iflag=0;
                    for (auto j=0;j<pdescenprogen[itimedescen][idescen].NumberofProgenitors;j++) {
                        if (pdescenprogen[itimedescen][idescen].dtoptype[j]==0) iflag=1;
                    }
                    if (iflag==1) continue;
                    //update by removing old links and adding new ones
                    RemoveLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                    dref[i]=dtemp[i];
                    AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                }
            }
        }
    }
    else if (opt.imultsteplinkcrit==MSLCMERITPRIMARYPROGEN) {
        for (Int_t i=0;i<numhalos;i++) {
            if (dref[i].NumberofDescendants==0 && dtemp[i].NumberofDescendants>0) {
                dref[i]=dtemp[i];
                AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
            }
            else if (dref[i].NumberofDescendants>0 && dtemp[i].NumberofDescendants>0) {
                if (dtemp[i].dtoptype[0]==0 && dref[i].dtoptype[0]!=0) {
                    //examine identified progenitor and see if it already has a 0 rank link
                    idescen=dtemp[i].DescendantList[0]-1;
                    itimedescen=dtemp[i].istep+itime;
                    iflag=0;
                    for (auto j=0;j<pdescenprogen[itimedescen][idescen].NumberofProgenitors;j++) {
                        if (pdescenprogen[itimedescen][idescen].dtoptype[j]==0) iflag=1;
                    }
                    if (iflag==1) continue;
                    RemoveLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                    dref[i]=dtemp[i];
                    AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                }
                else if (dtemp[i].dtoptype[0]<=dref[i].dtoptype[0] && dtemp[i].Merit[0]>dref[i].Merit[0]) {
                    RemoveLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                    dref[i]=dtemp[i];
                    AddLinksDescendantBasedProgenitorList(itime, i, dref[i], pdescenprogen);
                }
            }
        }
    }
}

///similar to construction of descendant list using the candidate progenitor list, but here working in reverse direction
///since objects can have many progenitors but an object should really have only one descendant, the ranking stored in progenindex
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
            pdescenprogen[did].progenindex.push_back(idescen);
            pdescenprogen[did].dtoptype.push_back(pdescen[k].dtoptype[idescen]);
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
    DescendantData *&pdescen, ProgenitorDataDescenBased *&pdescenprogen, int istep, Double_t meritlimit)
{
    PriorityQueue *pq;
    Int_t nattime,rank;
    Int_t progindex,descenindex,descenprogenindex;
    vector<int> progindexvec;
    vector<Double_t> meritvec;
    for (Int_t k=0;k<nhalos;k++) if (pdescenprogen[k].NumberofProgenitors>1)
    {
        //for each halo descandant rank its progenior haloes based on the merit at a given time
        nattime=0;
        for (auto iprogen=0;iprogen<pdescenprogen[k].NumberofProgenitors;iprogen++) nattime+=(pdescenprogen[k].deltat[iprogen]==istep);
        if (nattime<1) continue;
        //fill the priority queue so that the best merit are first
        pq=new PriorityQueue(nattime);
        for (auto iprogen=0;iprogen<pdescenprogen[k].NumberofProgenitors;iprogen++) if (pdescenprogen[k].deltat[iprogen]==istep)
        {
            pq->Push(iprogen,pdescenprogen[k].Merit[iprogen]);
        }
        //if the merit is poor then ranking starts at 1, otherwise, have a valid zero rank descendant/progenitor
        rank=(pq->TopPriority()<meritlimit);
        while(pq->Size()>0)
        {
            progindex=pq->TopQueue();
            descenindex=pdescenprogen[k].haloindex[progindex];
            descenprogenindex=pdescenprogen[k].progenindex[progindex];
            pdescenprogen[k].dtoptype[progindex]=rank;
            pdescen[descenindex].dtoptype[descenprogenindex]=rank;
            rank++;
            pq->Pop();
        }
        delete pq;
    }
}

void UpdateDescendantAcrossTimeUsingDescendantBasedProgenitorList(Int_t nhalos,
    DescendantData **&pdescenref, DescendantData *&pdescentemp, ProgenitorDataDescenBased *&pdescenprogen, int istep, Double_t meritlimit)
{
    PriorityQueue *pq;
    int iflag;
    Int_t rank;
    Int_t progindex,descenindex,descentemporalindex,descenprogenindex;
    vector<int> progindexvec;
    vector<Double_t> meritvec;
    for (Int_t k=0;k<nhalos;k++) if (pdescenprogen[k].NumberofProgenitors>1)
    {
        //if all objects are at the same time do nothing as this has been set already
        iflag=1;
        for (auto iprogen=1;iprogen<pdescenprogen[k].NumberofProgenitors;iprogen++) if (pdescenprogen[k].deltat[iprogen]!=pdescenprogen[k].deltat[0]) iflag=0;
        if (iflag) continue;
        //fill the priority queue so that the best merit are first
        pq=new PriorityQueue(pdescenprogen[k].NumberofProgenitors);
        for (auto iprogen=0;iprogen<pdescenprogen[k].NumberofProgenitors;iprogen++)
        {
            pq->Push(iprogen,pdescenprogen[k].Merit[iprogen]);
        }
        //if best rank is poor, lets start with 0
        rank=(pq->TopPriority()<meritlimit);
        while(pq->Size()>0)
        {
            progindex=pq->TopQueue();
            descenindex=pdescenprogen[k].haloindex[progindex];
            descentemporalindex=pdescenprogen[k].halotemporalindex[progindex];
            descenprogenindex=pdescenprogen[k].progenindex[progindex];
            pdescenprogen[k].dtoptype[progindex]=rank;
            //if deltat is that of the temporary descendant list, update the rank there
            if (pdescenprogen[k].deltat[progindex]==istep) pdescentemp[descenindex].dtoptype[descenprogenindex]=rank;
            //otherwise update rank of reference list
            else pdescenref[descentemporalindex][descenindex].dtoptype[descenprogenindex]=rank;
            rank++;
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
#ifdef USEMPI
        //for mpi communication need to store if anything is removed
        pdescenprogen[itimeprogen][did].removalhaloindex.push_back(pdescenprogen[itimeprogen][did].haloindex[k]);
        pdescenprogen[itimeprogen][did].removalhalotemporalindex.push_back(pdescenprogen[itimeprogen][did].halotemporalindex[k]);
#endif
        pdescenprogen[itimeprogen][did].haloindex.erase(pdescenprogen[itimeprogen][did].haloindex.begin()+k);
        pdescenprogen[itimeprogen][did].halotemporalindex.erase(pdescenprogen[itimeprogen][did].halotemporalindex.begin()+k);
        pdescenprogen[itimeprogen][did].Merit.erase(pdescenprogen[itimeprogen][did].Merit.begin()+k);
        pdescenprogen[itimeprogen][did].deltat.erase(pdescenprogen[itimeprogen][did].deltat.begin()+k);
        pdescenprogen[itimeprogen][did].progenindex.erase(pdescenprogen[itimeprogen][did].progenindex.begin()+k);
        pdescenprogen[itimeprogen][did].dtoptype.erase(pdescenprogen[itimeprogen][did].dtoptype.begin()+k);
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
        pdescenprogen[itimeprogen][did].progenindex.push_back(idescen);
        pdescenprogen[itimeprogen][did].dtoptype.push_back(pdescen.dtoptype[idescen]);
#ifdef USEMPI
        pdescenprogen[itimeprogen][did].MPITask.push_back(ThisTask);
#endif
        pdescenprogen[itimeprogen][did].NumberofProgenitors++;
    }
}

///Rank progenitors of a descandant based on their temporally generalized merit as viewed by their descendants.
void RankDescendantProgenitors(Int_t itimeprogen, HaloTreeData *&pht, ProgenitorDataDescenBased **&pdescenprogen, DescendantData **&pdescen, int iopttemporalmerittype)
{
    //if first snapshot can't have any progenitors
    if (itimeprogen==0) return;
    //if no halos before can't have any progenitors
    if (pht[itimeprogen-1].numhalos==0) return;
    PriorityQueue *pq;
    Int_t progindex,descenindex,descentemporalindex,descenprogenindex;
    Double_t generalizedmerit;
    unsigned short rank, maxrank;
    for (Int_t k=0;k<pht[itimeprogen].numhalos;k++) if (pdescenprogen[itimeprogen][k].NumberofProgenitors>1)
    {
        pq=new PriorityQueue(pdescenprogen[itimeprogen][k].NumberofProgenitors);
        /*
        maxrank=0;
        for (auto iprogen=0;iprogen<pdescenprogen[itimeprogen][k].NumberofProgenitors;iprogen++) {
            /*
            descenindex=pdescenprogen[itimeprogen][k].haloindex[iprogen];
            descentemporalindex=pdescenprogen[itimeprogen][k].halotemporalindex[iprogen];
            descenprogenindex=pdescenprogen[itimeprogen][k].progenindex[iprogen];
            rank=pdescen[descentemporalindex][descenindex].dtoptype[descenprogenindex];
            * /
            rank=pdescenprogen[itimeprogen][k].dtoptype[iprogen];
            if (maxrank<rank) maxrank=rank;
        }
        */
        for (auto iprogen=0;iprogen<pdescenprogen[itimeprogen][k].NumberofProgenitors;iprogen++)
        {
            //rank accoring to initial ranking and generalized merit.
/*
            descenindex=pdescenprogen[itimeprogen][k].haloindex[iprogen];
            descentemporalindex=pdescenprogen[itimeprogen][k].halotemporalindex[iprogen];
            descenprogenindex=pdescenprogen[itimeprogen][k].progenindex[iprogen];
            rank=pdescen[descentemporalindex][descenindex].dtoptype[descenprogenindex];
*/
            rank=pdescenprogen[itimeprogen][k].dtoptype[iprogen];
            generalizedmerit=pdescenprogen[itimeprogen][k].Merit[iprogen]/pow((Double_t)pdescenprogen[itimeprogen][k].deltat[iprogen],ALPHADELTAT);
            //the generalized merit takes rank into account
            //generalizedmerit/=(rank+1.0);
            generalizedmerit/=pow(4.0,rank);
            pq->Push(iprogen,generalizedmerit);
            //by using maxrank-rank+generalizedmerit, will have values staggared according to initial time specific rankings
            //pq->Push(iprogen,maxrank-rank+generalizedmerit);
        }
        rank=0;
        while(pq->Size()>0)
        {
            progindex=pq->TopQueue();
#ifdef USEMPI
            if (pdescenprogen[itimeprogen][k].MPITask[progindex]==ThisTask){
#endif
            descenindex=pdescenprogen[itimeprogen][k].haloindex[progindex];
            descentemporalindex=pdescenprogen[itimeprogen][k].halotemporalindex[progindex];
            descenprogenindex=pdescenprogen[itimeprogen][k].progenindex[progindex];
            pdescenprogen[itimeprogen][k].dtoptype[progindex]=rank;
            pdescen[descentemporalindex][descenindex].dtoptype[descenprogenindex]=rank;
#ifdef USEMPI
            }
#endif
            rank++;
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
    int mindtop,index,rank;
    Double_t maxmerit,generalizedmerit;
    unsigned long descen;
    vector<pair<Double_t,unsigned long>> merit;
    vector<int> dtop;
    vector<unsigned long> dindex;
    for (auto i=StartSnap;i<EndSnap-1;i++) {
        for (Int_t j=0;j<pht[i].numhalos;j++){
            if (pdescen[i][j].NumberofDescendants>1) {
                for (auto k=0;k<pdescen[i][j].NumberofDescendants;k++)
                {
                 	//rank accoring to initial ranking and generalized merit.
                    rank=pdescen[i][j].dtoptype[k];
                    generalizedmerit=pdescen[i][j].Merit[k];
                    generalizedmerit/=(rank+1.0);
                    //for ascending order merit, store 1/merit
                    merit.push_back(make_pair(1.0/generalizedmerit,k));
                    dtop.push_back(pdescen[i][j].dtoptype[k]);
                    dindex.push_back(pdescen[i][j].DescendantList[k]);
                }
                sort(merit.begin(),merit.end());
                for (auto k=0;k<pdescen[i][j].NumberofDescendants;k++) {
                    index=merit[k].second;
                  	pdescen[i][j].Merit[k]=1.0/merit[k].first;
                    pdescen[i][j].DescendantList[k]=dindex[index];
                    pdescen[i][j].dtoptype[k]=dtop[index];
                }
                merit.clear();
                dtop.clear();
                dindex.clear();
            }
        }
    }
}
//@}
//@}
