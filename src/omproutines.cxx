/*! \file omproutines.cxx
 *  \brief this file contains routines used with OpenMP compilation.

    OpenMP routines pertaining to complex OpenMP calls.
 */

#ifdef USEOPENMP

//-- For MPI

#include "stf.h"

/// \name routines which check to see if some search region overlaps with local mpi domain
//@{
///construct omp domains
OMP_Domain *OpenMPBuildDomains(Options &opt, const Int_t numompregions, KDTree *&tree, const Double_t rdist)
{
    OMP_Domain *ompdomain = new OMP_Domain[numompregions];
    Int_t noffset=0;
    Node *np;
    for (auto i=0;i<numompregions;i++) {
        np = (tree->FindLeafNode(noffset));
        ompdomain[i].ncount = np->GetCount();
        ompdomain[i].noffset = noffset;
        for (int j=0;j<3;j++) {
            ompdomain[i].bnd[j][0]=np->GetBoundary(j,0);
            ompdomain[i].bnd[j][1]=np->GetBoundary(j,1);
        }
        noffset+=ompdomain[i].ncount;
    }
    //determine for each omp region, what are neighbour omp regions within rdist
    Double_t xsearch[3][2];
    for (auto i=0;i<numompregions;i++) {
        for (auto k=0;k<3;k++) {xsearch[k][0]=ompdomain[i].bnd[k][0]-rdist;xsearch[k][1]=ompdomain[i].bnd[k][1]+rdist;}
        for (auto j=0;j<numompregions;j++) if (j!=i){
            if (OpenMPSearchForOverlap(xsearch,ompdomain[j].bnd, opt.p)) ompdomain[i].neighbour.push_back(j);
        }
        noffset+=ompdomain[i].ncount;
    }
    np=NULL;
    return ompdomain;
}

KDTree **OpenMPBuildLocalTrees(Options &opt, const Int_t numompregions, vector<Particle> &Part, OMP_Domain *ompdomain, Double_t *period)
{
    KDTree **tree3dfofomp = new KDTree*[numompregions];
    Int_t i;
    //get fof in each region
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(static) nowait
    for (i=0;i<numompregions;i++) {
        tree3dfofomp[i] = new KDTree(&Part.data()[ompdomain[i].noffset],ompdomain[i].ncount,opt.Bsize,tree3dfofomp[i]->TPHYS,tree3dfofomp[i]->KEPAN,100,0,0,0,period,NULL);
        tree3dfofomp[i]->OverWriteInputOrder();
    }
    }
    return tree3dfofomp;
}

///search if some region is in the local mpi domain
int OpenMPSearchForOverlap(Double_t xsearch[3][2], Double_t bnd[3][2], Double_t period){
    Double_t xsearchp[3][2];
    if (!((bnd[0][1] < xsearch[0][0]) || (bnd[0][0] > xsearch[0][1]) ||
        (bnd[1][1] < xsearch[1][0]) || (bnd[1][0] > xsearch[1][1]) ||
        (bnd[2][1] < xsearch[2][0]) || (bnd[2][0] > xsearch[2][1])))
            return 1;
    else {
        if (period==0) return 0;
        else {
            for (int j=0;j<3;j++) {xsearchp[j][0]=xsearch[j][0];xsearchp[j][1]=xsearch[j][1];}
            for (int j=0;j<3;j++) {
                if (!((bnd[j][1] < xsearch[j][0]+period) || (bnd[j][0] > xsearch[j][1]+period))) {xsearchp[j][0]+=period;xsearchp[j][1]+=period;}
                else if (!((bnd[j][1] < xsearch[j][0]-period) || (bnd[j][0] > xsearch[j][1]-period))) {xsearchp[j][0]-=period;xsearchp[j][1]-=period;}
            }
            if (!((bnd[0][1] < xsearchp[0][0]) || (bnd[0][0] > xsearchp[0][1]) ||
            (bnd[1][1] < xsearchp[1][0]) || (bnd[1][0] > xsearchp[1][1]) ||
            (bnd[2][1] < xsearchp[2][0]) || (bnd[2][0] > xsearchp[2][1])))
                return 1;
            else return 0;
        }
    }
}

/// Determine if a particle needs to be exported to another mpi domain based on a physical search radius
int OpenMPSearchForOverlap(Particle &Part, Double_t bnd[3][2], Double_t rdist, Double_t period)
{
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=Part.GetPosition(k)-rdist;xsearch[k][1]=Part.GetPosition(k)+rdist;}
    return OpenMPSearchForOverlap(xsearch,bnd,period);
}
///Determine if region fully in domain
int OpenMPInDomain(Double_t xsearch[3][2], Double_t bnd[3][2])
{
    return ((xsearch[0][0] > bnd[0][0]) && (xsearch[0][1] < bnd[0][1]) &&
    (xsearch[1][0] > bnd[1][0]) && (xsearch[1][1] < bnd[1][1]) &&
    (xsearch[2][0] > bnd[2][0]) && (xsearch[2][1] < bnd[2][1]));
}
///Determine if particle search region fully in domain
int OpenMPInDomain(Particle &Part, Double_t bnd[3][2], Double_t rdist)
{
  Double_t xsearch[3][2];
  for (auto k=0;k<3;k++) {xsearch[k][0]=Part.GetPosition(k)-rdist;xsearch[k][1]=Part.GetPosition(k)+rdist;}
  return OpenMPInDomain(xsearch,bnd);
}

Int_t OpenMPLocalSearch(Options &opt,
    const Int_t nbodies, vector<Particle> &Part, Int_t * &pfof, Int_t *&storeorgIndex,
    Int_tree_t *&Head, Int_tree_t *&Next,
    KDTree **&tree3dfofomp, Double_t *param, const Double_t rdist, const Int_t ompminsize, FOFcompfunc fofcmp,
    const Int_t numompregions, OMP_Domain *&ompdomain)
{
    Int_t i, orgIndex, ng = 0, ngtot=0;
    Int_t *p3dfofomp;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    auto time1=MyGetTime();
    cout<<ThisTask<<": Starting local openmp searches "<<endl;
    #pragma omp parallel default(shared) \
    private(i,p3dfofomp,orgIndex, ng)
    {
    #pragma omp for schedule(dynamic) nowait reduction(+:ngtot)
    for (i=0;i<numompregions;i++) {
        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1) p3dfofomp=tree3dfofomp[i]->FOFCriterionSetBasisForLinks(fofcmp,param,ng,ompminsize,0,0,FOFchecktype, &Head[ompdomain[i].noffset], &Next[ompdomain[i].noffset]);
        else p3dfofomp=tree3dfofomp[i]->FOF(rdist,ng,ompminsize,0, &Head[ompdomain[i].noffset], &Next[ompdomain[i].noffset]);
        if (ng > 0) {
            for (int j=ompdomain[i].noffset;j<ompdomain[i].noffset+ompdomain[i].ncount;j++)
            if (p3dfofomp[Part[j].GetID()]>0)
            {
                orgIndex = storeorgIndex[Part[j].GetID()+ompdomain[i].noffset];
                pfof[orgIndex] = p3dfofomp[Part[j].GetID()]+ompdomain[i].noffset;
            }
        }
        delete[] p3dfofomp;
        ompdomain[i].numgroups = ng;
        ngtot += ng;
    }
    }
    cout<<ThisTask<<" finished local search "<<ngtot<<" in "<<MyElapsedTime(time1)<<endl;
    return ngtot;
}

///Saerch particles to see if they overlap other OpenMP domains
OMP_ImportInfo *OpenMPImportParticles(Options &opt, const Int_t nbodies, vector<Particle> &Part,
    Int_t * &pfof, Int_t *&storeorgIndex,
    const Int_t numompregions, OMP_Domain *&ompdomain, const Double_t rdist,
    Int_t *&omp_nrecv_total, Int_t *&omp_nrecv_offset, Int_t &importtotal)
{
    Int_t i,j,orgIndex;
    importtotal=0;
    auto time1=MyGetTime();
    OMP_ImportInfo *ompimport;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    cout<<ThisTask<<": Starting import build "<<endl;
    for (i=0;i<numompregions;i++) omp_nrecv_total[i]=omp_nrecv_offset[i]=0;
    for (i=0;i<numompregions;i++) {
        for (j=ompdomain[i].noffset;j<ompdomain[i].noffset+ompdomain[i].ncount;j++) {
            if (OpenMPInDomain(Part[j],ompdomain[i].bnd,rdist)) continue;
            for (auto k: ompdomain[i].neighbour) {
                if (OpenMPSearchForOverlap(Part[j],ompdomain[k].bnd,rdist,opt.p))  omp_nrecv_total[k] += 1;
            }
        }
    }

    for (i=1;i<numompregions;i++) omp_nrecv_offset[i] = omp_nrecv_offset[i-1]+omp_nrecv_total[i-1];
    for (i=0;i<numompregions;i++) {
        if (opt.iverbose > 1) cout<<ThisTask<<" omp region "<<i<<" is importing "<<omp_nrecv_total[i]<<endl;
        importtotal += omp_nrecv_total[i];
        omp_nrecv_total[i] = 0;
    }
    if (importtotal == 0) {ompimport=NULL; return ompimport;}

    ompimport = new OMP_ImportInfo[importtotal];
    for (i=0;i<numompregions;i++) {
        for (j=ompdomain[i].noffset;j<ompdomain[i].noffset+ompdomain[i].ncount;j++) {
            if (OpenMPInDomain(Part[j],ompdomain[i].bnd,rdist)) continue;
            for (auto k: ompdomain[i].neighbour) {
                if (OpenMPSearchForOverlap(Part[j],ompdomain[k].bnd,rdist,opt.p)) {
                    orgIndex = storeorgIndex[Part[j].GetID()+ompdomain[i].noffset];
                    ompimport[omp_nrecv_total[k]+omp_nrecv_offset[k]].index = j;
                    ompimport[omp_nrecv_total[k]+omp_nrecv_offset[k]].pfof = pfof[orgIndex];
                    ompimport[omp_nrecv_total[k]+omp_nrecv_offset[k]].task = i;
                    omp_nrecv_total[k] += 1;
                }
            }
        }
    }
    cout<<ThisTask<<" finished import "<<MyElapsedTime(time1)<<endl;
    return ompimport;
}

void OpenMPLinkAcross(Options &opt,
    Int_t nbodies, vector<Particle> &Part, Int_t * &pfof, Int_t *&storeorgIndex,
    Int_tree_t *&Head, Int_tree_t *&Next,
    Double_t *param, FOFcheckfunc &fofcheck,
    const Int_t numompregions, OMP_Domain *&ompdomain, KDTree **tree3dfofomp,
    Int_t *&omp_nrecv_total, Int_t *&omp_nrecv_offset, OMP_ImportInfo* &ompimport)
{
    //now begin linking across
    Int_t i;
    Int_t omp_links_across_total, numloops;
    Int_t nt, orgIndex, curIndex, *nn=new Int_t[nbodies], pfofcomp;
    Int_t maxgroupnum;
    vector<Int_t> localnewgroup(numompregions,0);
    int curTask;
    Coordinate x;
    auto time1=MyGetTime();
    Particle *Pval;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    if (NProcs > 1) {
        maxgroupnum = ompdomain[numompregions-1].numgroups + ompdomain[numompregions-1].noffset;
    }

    cout<<ThisTask<<": Starting linking across OpenMP domains"<<endl;
    numloops = 0;
    do {
        omp_links_across_total = 0;
        #pragma omp parallel default(shared) \
        private(i, orgIndex, curIndex, x, nt, Pval, pfofcomp)
        {
        #pragma omp for nowait reduction(+:omp_links_across_total)
        for (i=0;i<numompregions;i++) {
            for (auto j=0;j<omp_nrecv_total[i];j++) {
                Pval=&Part[ompimport[omp_nrecv_offset[i]+j].index];
                pfofcomp = ompimport[omp_nrecv_offset[i]+j].pfof;
                /// before ignored pfofcomp but since cannot guarantee openmp
                /// will be complete if MPI is used.
                if (pfofcomp == 0 && NProcs == 1) continue;
                //for each imported particle, find all particles within search window
                for (auto k=0;k<3;k++) x[k]=Pval->GetPosition(k);
                nt=tree3dfofomp[i]->SearchBallPosTagged(x, param[1], &nn[ompdomain[i].noffset]);
                for (auto k=0;k<nt;k++) {
                    curIndex=nn[k+ompdomain[i].noffset]+ompdomain[i].noffset;
                    //check that at least on of the particles meets the type criterion if necessary
                    if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                        if (fofcheck(Part[curIndex],param)!=0 && fofcheck(*Pval,param)!=0) continue;

                    orgIndex = storeorgIndex[Part[curIndex].GetID()+ompdomain[i].noffset];
                    //otherwise, change these particles to local group id if local group id smaller
                    //if local particle in a group
                    if (pfof[orgIndex]>0 && pfofcomp > 0)  {
                        //only change if both particles are appropriate type and group ids indicate local needs to be exported
                        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                            if (!(fofcheck(Part[curIndex],param)==0 && fofcheck(*Pval,param)==0)) continue;
                        //if local group id is larger, change locally
                        if (pfof[orgIndex] > pfofcomp) {
                            Int_t ss = Head[nn[k+ompdomain[i].noffset]+ompdomain[i].noffset];
                            do{
                                orgIndex = storeorgIndex[Part[ss+ompdomain[i].noffset].GetID()+ompdomain[i].noffset];
                                pfof[orgIndex]=pfofcomp;
                            }while((ss = Next[ss+ompdomain[i].noffset]) >= 0);
                            omp_links_across_total++;
                        }
                    }
                    //if local particle not in a group and export is appropriate type, link
                    else {
                        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                            if (fofcheck(*Pval,param)!=0) continue;
                        //if local particle not in group and associated particle in group
                        //set the local particle to the appropriate id
                        if (pfofcomp > 0) {
                            pfof[orgIndex] = pfofcomp;
                        }
                        ///otherwise, if at this point, omp domains are incomplete
                        ///as volume also split in mpi, so increase number of local
                        ///groups and assign particle to this new group 
                        else {
                            localnewgroup[i]++;
                            pfof[orgIndex] = maxgroupnum + ompdomain[i].noffset + localnewgroup[i];
                        }
                        omp_links_across_total++;
                    }
                }
            }
        }
        }

        #pragma omp parallel default(shared) \
        private(i, orgIndex, curIndex, curTask)
        {
        #pragma omp for nowait
        for (i=0;i<numompregions;i++) {
            for (auto j=0;j<omp_nrecv_total[i];j++) {
                curIndex=ompimport[omp_nrecv_offset[i]+j].index;
                curTask=ompimport[omp_nrecv_offset[i]+j].task;
                orgIndex=storeorgIndex[Part[curIndex].GetID()+ompdomain[curTask].noffset];
                ompimport[omp_nrecv_offset[i]+j].pfof=pfof[orgIndex];
            }
        }
        }
        numloops++;
        if (opt.iverbose>1) cout<<ThisTask<<" linking across at loop "<<numloops<<" having found "<<omp_links_across_total<<" links"<<endl;
    }while(omp_links_across_total>0);
    delete[] nn;
    cout<<ThisTask<<" finished linking "<<MyElapsedTime(time1)<<endl;
}

Int_t OpenMPResortParticleandGroups(Int_t nbodies, vector<Particle> &Part, Int_t *&pfof, Int_t minsize)
{
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    Int_t ngroups = 0, newnumgroups = 0;;
    vector<Int_t> numingroup, index;
    unordered_map<Int_t, Int_t> pfofoldtonew;
    PriorityQueue *pq;

    //init data
    for (auto i=0;i<nbodies;i++) if (ngroups < pfof[i]) ngroups = pfof[i];
    numingroup.resize(ngroups+1,0);
    index.resize(ngroups+1);
    for (auto i=0;i<=ngroups;i++) index[i] = i;
    for (auto i=0;i<nbodies;i++) if (pfof[i]>0) numingroup[pfof[i]]++;
    for (auto i=1;i<=ngroups;i++) {
        if (numingroup[i]<minsize) numingroup[i] = 0;
        else newnumgroups++;
    }
    //if no groups are large enough, zero and return
    if (newnumgroups == 0) {
        #pragma omp parallel for default(shared)
        for (auto i=0;i<nbodies;i++) {
            pfof[i] = 0;
        }
        ngroups = newnumgroups;
        return ngroups;
    }

    //otherwise, remap group ids so as to be in decreasing group size
    //first set all pfof values for all groups below minsize to zero
    for (auto i=0;i<nbodies;i++) {
        if (pfof[i] == 0) continue;
        if (numingroup[pfof[i]]<minsize) pfof[i] = 0;
    }
    pq = new PriorityQueue(newnumgroups);
    for (auto i=1;i<=ngroups;i++) {
        if (numingroup[i] == 0) continue;
        pq->Push(index[i],numingroup[i]);
    }
    ngroups = newnumgroups;
    newnumgroups = 0;
    //generate map
    for (auto i=0;i<ngroups;i++) {
        pfofoldtonew[pq->TopQueue()] = ++newnumgroups;
        pq->Pop();
    }
    delete pq;
    //set new group id values stored in pfof
    for (auto i=0;i<nbodies;i++) {
        if (pfof[i] == 0) continue;
        pfof[i]=pfofoldtonew[pfof[i]];
    }
    return ngroups;
}

void OpenMPHeadNextUpdate(const Int_t nbodies, vector<Particle> &Part, const Int_t numgroups, Int_t *&pfof, Int_tree_t *&Head, Int_tree_t *&Next){
    Int_t *numingroup;
    Int_t **pglist;
    numingroup=BuildNumInGroup(nbodies, numgroups, pfof);
    pglist=BuildPGList(nbodies, numgroups, numingroup, pfof, Part.data());
    for (auto i=0;i<nbodies;i++) {
        Head[i]=i;
        Next[i]=-1;
    }
    for (auto i=1;i<=numgroups;i++) {
        for (auto j=1;j<numingroup[i];j++) Head[pglist[i][j]]=pglist[i][0];
        for (auto j=0;j<numingroup[i]-1;j++) Next[pglist[i][j]]=pglist[i][j+1];
    }
    for (auto i=1;i<=numgroups;i++) delete[] pglist[i];
    delete[] numingroup;
    delete[] pglist;
}

//@}

#endif
