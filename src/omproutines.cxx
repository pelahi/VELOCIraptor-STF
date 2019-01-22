/*! \file omproutines.cxx
 *  \brief this file contains routines used with OpenMP compilation.

    OpenMP routines pertaining to complex OpenMP calls.
 */

#ifdef USEOPENMP

//-- For MPI

#include "stf.h"

/// \name routines which check to see if some search region overlaps with local mpi domain
//@{
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

void OpenMPLinkAcross(Options &opt,
    Int_t nbodies, vector<Particle> &Part, Int_t * &pfof,
    Int_t *&storetype, Int_tree_t *&Head, Int_tree_t *&Next,
    Double_t *param, FOFcheckfunc &fofcheck,
    const Int_t numompregions, OMP_Domain *&ompdomain, KDTree **tree3dfofomp,
    Int_t *&omp_nrecv_total, Particle** &Partompimport)
{
    //now begin linking across
    Int_t i;
    Int_t omp_links_across_total;
    Int_t nt, orgIndex, curIndex, *nn=new Int_t[nbodies];
    int numloops = 0;
    Double_t rdist=sqrt(param[1]), rdist2=param[1];
    Coordinate x;
#ifndef USEMPI
    int ThisTask =0;
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
                for (auto k=0;k<3;k++) x[k]=Partompimport[i][j].GetPosition(k);
                nt=tree3dfofomp[i]->SearchBallPosTagged(x, rdist2, &nn[ompdomain[i].noffset]);
                for (auto k=0;k<nt;k++) {
                    curIndex=nn[k+ompdomain[i].noffset]+ompdomain[i].noffset;
                    //check that at least on of the particles meets the type criterion if necessary
                    if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                        if (fofcheck(Part[curIndex],param)!=0 && fofcheck(Partompimport[i][j],param)!=0) continue;
                    orgIndex = storetype[Part[curIndex].GetID()+ompdomain[i].noffset];
                    //otherwise, change these particles to local group id if local group id smaller
                    //if local particle in a group
                    if (pfof[orgIndex]>0)  {
                        //only change if both particles are appropriate type and group ids indicate local needs to be exported
                        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                            if (!(fofcheck(Part[curIndex],param)==0 && fofcheck(Partompimport[i][j],param)==0)) continue;
                        //if local group id is larger, change locally
                        if(pfof[orgIndex] > Partompimport[i][j].GetID()) {
                            Int_t ss = Head[nn[k+ompdomain[i].noffset]+ompdomain[i].noffset];
                            do{
                                orgIndex = storetype[Part[ss+ompdomain[i].noffset].GetID()+ompdomain[i].noffset];
                                pfof[orgIndex]=Partompimport[i][j].GetID();
                            }while((ss = Next[ss+ompdomain[i].noffset]) >= 0);
                            omp_links_across_total++;
                        }
                    }
                    //if local particle not in a group and export is appropriate type, link
                    else {
                        if (opt.partsearchtype==PSTALL && opt.iBaryonSearch>1)
                            if (fofcheck(Partompimport[i][j],param)!=0) continue;
                        pfof[orgIndex]=Partompimport[i][j].GetID();
                        omp_links_across_total++;
                    }
                }

            }
        }
        }
        numloops++;
        if (opt.iverbose > 1) cout<<ThisTask<<" at loop "<<numloops<<" having found "<<omp_links_across_total<<" links across omp domains"<<endl;
    }while(omp_links_across_total>0);
    delete[] nn;
}

Int_t OpenMPResortParticleandGroups(Int_t nbodies, vector<Particle> &Part, Int_t *&pfof, Int_t minsize)
{
    Int_t start, ngroups=0;
    Int_t *numingroup, **plist;

    //now get number of groups and reorder group ids
    Int_t *storepid = new Int_t[nbodies];
    for (auto i=0;i<nbodies;i++) {
        storepid[i] = Part[i].GetPID();
        Part[i].SetPID(-pfof[i]);
    }
    qsort(Part.data(),nbodies,sizeof(Particle),PIDCompare);

    //determine the # of groups, their size and the current group ID
    for (auto i=0,start=0;i<nbodies;i++) {
        if (Part[i].GetPID()!=Part[start].GetPID()) {
            //if group is too small set pid to zero, which currently is used to store the group id
            if ((i-start)<minsize) for (Int_t j=start;j<i;j++) Part[j].SetPID(0);
            else ngroups++;
            start=i;
        }
        if (Part[i].GetPID()==0) break;
    }

    //again resort to move untagged particles to the end.
    qsort(Part.data(),nbodies,sizeof(Particle),PIDCompare);

    //now adjust pfof and ids.
    for (auto i=0;i<nbodies;i++) {pfof[Part[i].GetID()]=-Part[i].GetPID();Part[i].SetPID(storepid[Part[i].GetID()]);}
    delete[] storepid;

    numingroup=new Int_t[ngroups+1];
    plist=new Int_t*[ngroups+1];
    ngroups=1;//offset as group zero is untagged
    for (auto i=0,start=0;i<nbodies;i++) {
        if (pfof[Part[i].GetID()]!=pfof[Part[start].GetID()]) {
            numingroup[ngroups]=i-start;
            plist[ngroups]=new Int_t[numingroup[ngroups]];
            for (auto j=start,count=0;j<i;j++) plist[ngroups][count++]=Part[j].GetID();
            ngroups++;
            start=i;
        }
        if (pfof[Part[i].GetID()]==0) break;
    }
    ngroups--;

    //reorder groups ids according to size
    qsort(Part.data(),nbodies,sizeof(Particle), IDCompare);
    ReorderGroupIDs(ngroups,ngroups,numingroup,pfof,plist);
    for (auto i=1;i<=ngroups;i++) delete[] plist[i];
    delete[] plist;
    delete[] numingroup;
    return ngroups;
}



//@}

#endif
