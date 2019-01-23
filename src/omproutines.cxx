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
            ompdomain[i].bnd[j][0]=np->GetBoundary(j,0)*0.99;
            ompdomain[i].bnd[j][1]=np->GetBoundary(j,1)*1.01;
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
    #pragma omp for schedule(dynamic) nowait
    for (i=0;i<numompregions;i++) {
        tree3dfofomp[i] = new KDTree(&Part.data()[ompdomain[i].noffset],ompdomain[i].ncount,opt.Bsize,tree3dfofomp[i]->TPHYS,tree3dfofomp[i]->KEPAN,100,0,0,0,period);
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
int OpenMPSearchForOverlap(Particle &Part, Double_t bnd[3][2], Double_t rdist, Double_t period){
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=Part.GetPosition(k)-rdist;xsearch[k][1]=Part.GetPosition(k)+rdist;}
    return OpenMPSearchForOverlap(xsearch,bnd,period);
}

///Saerch particles to see if they overlap other OpenMP domains
Particle *OpenMPImportParticles(Options &opt, const Int_t nbodies, vector<Particle> &Part, Int_t * &pfof, Int_t *&storetype,
    const Int_t numompregions, OMP_Domain *&ompdomain, const Double_t rdist,
    Int_t *&omp_nrecv_total, Int_t *&omp_nrecv_offset)
{
    Int_t i,orgIndex;
    int omptask;
    Int_t importtotal=0;
    double time1=MyGetTime();
    Particle *Partompimport;
    Int_t *Partompimportindex;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    cout<<ThisTask<<": Starting import build "<<endl;
    for (i=0;i<numompregions;i++) omp_nrecv_total[i]=omp_nrecv_offset[i]=0;
    #pragma omp parallel default(shared) \
    private(i,orgIndex,omptask)
    {
    #pragma omp for
    for (i=0;i<numompregions;i++) {
        for (auto j=ompdomain[i].noffset;j<ompdomain[i].noffset+ompdomain[i].ncount;j++) {
            orgIndex = storetype[Part[j].GetID()+ompdomain[i].noffset];
            if (pfof[orgIndex] == 0) continue;
            //for (auto k=0;k<numompregions;k++) if (k!=i) {
            for (auto k: ompdomain[i].neighbour) {
                if (OpenMPSearchForOverlap(Part[i],ompdomain[k].bnd,rdist,opt.p)) omp_nrecv_total[k] += 1;
            }
        }
    }
    }

    for (i=1;i<numompregions;i++) omp_nrecv_offset[i] = omp_nrecv_offset[i-1]+omp_nrecv_total[i-1];
    for (i=0;i<numompregions;i++) {
        if (opt.iverbose > 1) cout<<ThisTask<<" omp region "<<i<<" is importing "<<omp_nrecv_total[i]<<endl;
        importtotal += omp_nrecv_total[i];
        omp_nrecv_total[i] = 0;
    }
    Partompimport = new Particle[importtotal];

    #pragma omp parallel default(shared) \
    private(i,orgIndex,omptask)
    {
    #pragma omp for
    for (i=0;i<numompregions;i++) {
        for (auto j=ompdomain[i].noffset;j<ompdomain[i].noffset+ompdomain[i].ncount;j++) {
            orgIndex = storetype[Part[j].GetID()+ompdomain[i].noffset];
            if (pfof[orgIndex] == 0) continue;
            //for (auto k=0;k<numompregions;k++) if (k!=i) {
            for (auto k: ompdomain[i].neighbour) {
                if (OpenMPSearchForOverlap(Part[i],ompdomain[k].bnd,rdist,opt.p)) {
                    Partompimport[omp_nrecv_total[k]+omp_nrecv_offset[k]] = Part[j];
                    Partompimport[omp_nrecv_total[k]+omp_nrecv_offset[k]].SetID(pfof[orgIndex]);
                    omp_nrecv_total[k] += 1;
                }
            }
        }
    }
    }
    cout<<ThisTask<<" finished import "<<MyGetTime()-time1<<endl;
    return Partompimport;
}

void OpenMPLinkAcross(Options &opt,
    Int_t nbodies, vector<Particle> &Part, Int_t * &pfof,
    Int_t *&storetype, Int_tree_t *&Head, Int_tree_t *&Next,
    Double_t *param, FOFcheckfunc &fofcheck,
    const Int_t numompregions, OMP_Domain *&ompdomain, KDTree **tree3dfofomp,
    Int_t *&omp_nrecv_total, Int_t *&omp_nrecv_offset, Particle* &Partompimport)
{
    //now begin linking across
    Int_t i;
    Int_t omp_links_across_total;
    Int_t nt, orgIndex, curIndex, *nn=new Int_t[nbodies];
    Coordinate x;
    double time1=MyGetTime();
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    cout<<ThisTask<<": Starting linking across OpenMP domains"<<endl;
    do {
        omp_links_across_total = 0;
        #pragma omp parallel default(shared) \
        private(i,orgIndex,curIndex, x, nt)
        {
        #pragma omp for schedule(dynamic,1) nowait reduction(+:omp_links_across_total)
        for (i=0;i<numompregions;i++) {
            for (auto j=0;j<omp_nrecv_total[i];j++) {
                //for each imported particle, find all particles within search window
                for (auto k=0;k<3;k++) x[k]=Partompimport[omp_nrecv_offset[i]+j].GetPosition(k);
                nt=tree3dfofomp[i]->SearchBallPosTagged(x, param[1], &nn[ompdomain[i].noffset]);
                for (auto k=0;k<nt;k++) {
                    curIndex=nn[k+ompdomain[i].noffset]+ompdomain[i].noffset;
                    //check that at least on of the particles meets the type criterion if necessary
                    if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                        if (fofcheck(Part[curIndex],param)!=0 && fofcheck(Partompimport[omp_nrecv_offset[i]+j],param)!=0) continue;

                    orgIndex = storetype[Part[curIndex].GetID()+ompdomain[i].noffset];
                    //otherwise, change these particles to local group id if local group id smaller
                    //if local particle in a group
                    if (pfof[orgIndex]>0)  {
                        //only change if both particles are appropriate type and group ids indicate local needs to be exported
                        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                            if (!(fofcheck(Part[curIndex],param)==0 && fofcheck(Partompimport[omp_nrecv_offset[i]+j],param)==0)) continue;
                        //if local group id is larger, change locally
                        if(pfof[orgIndex] > Partompimport[omp_nrecv_offset[i]+j].GetPID()) {
                            Int_t ss = Head[nn[k+ompdomain[i].noffset]+ompdomain[i].noffset];
                            do{
                                orgIndex = storetype[Part[ss+ompdomain[i].noffset].GetID()+ompdomain[i].noffset];
                                pfof[orgIndex]=Partompimport[omp_nrecv_offset[i]+j].GetPID();
                            }while((ss = Next[ss+ompdomain[i].noffset]) >= 0);
                            omp_links_across_total++;
                        }
                    }
                    //if local particle not in a group and export is appropriate type, link
                    else {
                        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                            if (fofcheck(Partompimport[omp_nrecv_offset[i]+j],param)!=0) continue;
                        pfof[orgIndex]=Partompimport[omp_nrecv_offset[i]+j].GetPID();
                        omp_links_across_total++;
                    }
                }
            }
        }
        }
    }while(omp_links_across_total>0);
    delete[] nn;
    cout<<ThisTask<<" finished linking "<<MyGetTime()-time1<<endl;
}

Int_t OpenMPResortParticleandGroups(Int_t nbodies, vector<Particle> &Part, Int_t *&pfof, Int_t minsize)
{
    Int_t start, ngroups=0;
    Int_t *numingroup, **plist;
    //now get number of groups and reorder group ids
    for (auto i=0;i<nbodies;i++) Part[i].SetID(-pfof[i]);
    //used to use ID store store group id info
    qsort(Part.data(),nbodies,sizeof(Particle),IDCompare);

    //determine the # of groups, their size and the current group ID
    for (auto i=0,start=0;i<nbodies;i++) {
        if (Part[i].GetID()!=Part[start].GetID()) {
            //if group is too small set type to zero, which currently is used to store the group id
            if ((i-start)<minsize) for (Int_t j=start;j<i;j++) Part[j].SetID(0);
            else ngroups++;
            start=i;
        }
        if (Part[i].GetID()==0) break;
    }
    //again resort to move untagged particles to the end.
    qsort(Part.data(),nbodies,sizeof(Particle),IDCompare);

    //now adjust pfof and ids.
    for (auto i=0;i<nbodies;i++) {pfof[i]=-Part[i].GetID();Part[i].SetID(i);}
    numingroup=new Int_t[ngroups+1];
    plist=new Int_t*[ngroups+1];
    ngroups=1;//offset as group zero is untagged
    for (auto i=0,start=0;i<nbodies;i++) {
        if (pfof[i]!=pfof[start]) {
            numingroup[ngroups]=i-start;
            plist[ngroups]=new Int_t[numingroup[ngroups]];
            for (auto j=start,count=0;j<i;j++) plist[ngroups][count++]=j;
            ngroups++;
            start=i;
        }
        if (pfof[i]==0) break;
    }
    ngroups--;

    //reorder groups ids according to size
    ReorderGroupIDs(ngroups,ngroups,numingroup,pfof,plist);
    for (auto i=1;i<=ngroups;i++) delete[] plist[i];
    delete[] plist;
    delete[] numingroup;
    return ngroups;
}



//@}

#endif
