/*! \file unbind.cxx
 *  \brief this file contains routines to check if groups are self-bound and if not unbind them as requried

    \todo Need to improve the gravity calculation (ie: for tree potential use something other than just monopole and also apply corrections if necessary).
    \todo Need to clean up unbind proceedure, ensure its mpi compatible and can be combined with a pglist output easily
 */

#include "stf.h"
#include <random>

///\name Tree-Potential routines
//@{
///subroutine that generates node list for tree gravity calculation
void GetNodeList(Node *np, Int_t &ncell, Node **nodelist, const Int_t bsize){
    nodelist[ncell]=np;
    nodelist[ncell]->SetID(ncell);
    if (np->GetCount()>bsize){
        ncell++;GetNodeList(((SplitNode*)np)->GetLeft(),ncell,nodelist,bsize);
        ncell++;GetNodeList(((SplitNode*)np)->GetRight(),ncell,nodelist,bsize);
    }
    //else ncell++;
}

///subroutine that marks a cell for a given particle in tree-walk
inline void MarkCell(Node *np, Int_t *marktreecell, Int_t *markleafcell, Int_t &ntreecell, Int_t &nleafcell, Double_t *r2val, const Int_t bsize, Double_t *cR2max, Coordinate *cm, Double_t *cmtot, Coordinate xpos, Double_t eps2){
    Int_t nid=np->GetID();
    Double_t r2;
    r2=0;
    //determine the distance from cells cm to particle
    for (int k=0;k<3;k++)r2+=(cm[nid][k]-xpos[k])*(cm[nid][k]-xpos[k]);
    //if the particle is not distant enough to treat cell (and all subcells) as a mono (or quad) pole mass then
    //enter the cell so long as it is not a leaf node (or minimum cell size)
    //if it is a minimum cell size, the cell is marked with a ileafflag
    if (r2<cR2max[nid]) {
        if (np->GetCount()>bsize){
            MarkCell(((SplitNode*)np)->GetLeft(),marktreecell,markleafcell,ntreecell,nleafcell,r2val,bsize,cR2max,cm,cmtot,xpos,eps2);
            MarkCell(((SplitNode*)np)->GetRight(),marktreecell,markleafcell,ntreecell,nleafcell,r2val,bsize,cR2max,cm,cmtot,xpos,eps2);
        }
        else markleafcell[nleafcell++]=nid;
    }
    else {
        r2val[ntreecell]=cmtot[nid]/sqrt(r2+eps2);
        marktreecell[ntreecell++]=nid;
    }
}

//@}

//@{
inline bool CheckGroupForBoundness(Options &opt, Double_t &Efrac, Double_t &maxE, Int_t ning) {
    bool unbindcheck=false;
    if (opt.uinfo.unbindtype==USYSANDPART) {
        if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(ning>=opt.MinSize)) unbindcheck=true;
        else unbindcheck=false;
    }
    else if (opt.uinfo.unbindtype==UPART) {
        if ((maxE>0)&&(ning>=opt.MinSize))unbindcheck=true;
        else unbindcheck=false;
    }
    return unbindcheck;
}

inline void FillUnboundArrays(Options &opt, int maxunbindsize,
    Int_t ning, Particle *groupPart, const Double_t &Efrac,
    Int_t *&nEplusid, int *&Eplusflag, Int_t &nEplus, bool &unbindcheck
)
{
    nEplus=0;
    if (unbindcheck==false) return;
    //if just looking at particle then add to removal list till energy >0
    if (opt.uinfo.unbindtype==UPART) {
        for (auto j=0;j<maxunbindsize;j++) {
            if (groupPart[ning-1-j].GetDensity()>0) {
                nEplusid[nEplus++]=ning-1-j;
                Eplusflag[ning-1-j]=1;
            }
            else break;
        }
    }
    //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
    else if (opt.uinfo.unbindtype==USYSANDPART) {
        Int_t nEfrac=0;
        if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*ning;
        for (auto j=0;j<maxunbindsize;j++) {
            if (groupPart[ning-1-j].GetDensity()>0 || nEplus<nEfrac)
            {
                nEplusid[nEplus++]=ning-1-j;
                Eplusflag[ning-1-j]=1;
            }
            else break;
        }
    }
    //if number of unbound particles is much less than total number of particles
    //which can be unbound in one go, set unbinding check to false
    if (nEplus < opt.uinfo.maxallowedunboundfrac*ning) {
        nEplus = 0;
        unbindcheck = false;
    }
}

///remove particles deemed unbound
inline void RemoveUnboundParticles(Int_t i, Int_t *&pfof, Int_t &ning, Int_t *&pglist, Particle *groupPart,
            Int_t &nEplus, Int_t *&nEplusid, int *&Eplusflag)
{
    //adjust the pfof array and move unbound particles to the end of the array
    for (auto j=0;j<nEplus;j++) pfof[groupPart[nEplusid[j]].GetPID()]=0;
    for (auto j=0;j<nEplus;j++) groupPart[nEplusid[j]].SetPID(-1);
    ning-=nEplus;
}

inline void RemoveUnboundParticles(Int_t i, Int_t *&pfof, Int_t &ning, Particle *groupPart,
            Int_t &nEplus, Int_t *&nEplusid, int *&Eplusflag)
{
    //adjust the pfof array and move unbound particles to the end of the array
    for (auto j=0;j<nEplus;j++) pfof[groupPart[nEplusid[j]].GetID()]=0;
    ning-=nEplus;
}

///correct cm for removal of all least bound particles
inline void UpdateCMForUnboundParticles(Options &opt,
    Double_t &mass, Coordinate &cmvel,
    Int_t &nig, Particle *groupPart,
    Int_t &nEplus, Int_t *&nEplusid, int *&Eplusflag)
{
    double temp=1.0/mass, temp2=0.;
    if (opt.uinfo.fracpotref == 1.0) {
        for (auto j=0;j<nEplus;j++) {
            for (auto k=0;k<3;k++)
                cmvel[k]-=groupPart[nEplusid[j]].GetVelocity(k)*groupPart[nEplusid[j]].GetMass()*temp;
            temp2+=groupPart[nEplusid[j]].GetMass();
        }
        temp=mass/(mass-temp2);
        cmvel*=temp;
        mass-=temp2;
    }
}

/// Update the potential if necessary for small groups
inline void UpdatePotentialForUnboundParticlesPP(Options &opt,
    Int_t &nig, Particle *groupPart,
    Int_t &nEplus, Int_t *&nEplusid, int *&Eplusflag)
{
    //if keeping background then do nothing
    if (opt.uinfo.bgpot!=0) return;
    Double_t r2, poti, eps2=opt.uinfo.eps*opt.uinfo.eps;
#ifdef NOMASS
    Double_t mv2=opt.MassValue*opt.MassValue;
#endif
    for (auto k=0;k<nEplus;k++) {
        for (auto j=0;j<nig;j++) {
            if (j!=nEplusid[k]) {
                r2=0.;for (auto n=0;n<3;n++) r2+=pow(groupPart[nEplusid[k]].GetPosition(n)-groupPart[j].GetPosition(n),2.0);
                r2+=eps2;
                r2=1.0/sqrt(r2);
                poti=groupPart[j].GetPotential()+opt.G*(groupPart[nEplusid[k]].GetMass()*groupPart[j].GetMass())*r2;
#ifdef NOMASS
                poti*=mv2;
#endif
                groupPart[j].SetPotential(poti);
            }
        }
    }
}

/// Update the potential if necessary for large groups
inline void UpdatePotentialForUnboundParticles(Options &opt,
    Int_t &nig, Particle *groupPart,
    Int_t &nEplus, Int_t *&nEplusid, int *&Eplusflag)
{
    int iunbindsizeflag;
    Double_t r2, poti, eps2=opt.uinfo.eps*opt.uinfo.eps;
#ifdef NOMASS
    Double_t mv2=opt.MassValue*opt.MassValue;
#endif

    if (opt.uinfo.bgpot!=0) return;
    //if ignore the background then adjust the potential energy of the particles
    //for large groups with many particles removed more computationally effective to simply
    //recalculate the potential energy after removing particles
    //for smaller number of particles removed, simply remove the contribution of this particle
    //from all others. The change in efficiency occurs at roughly nEplus>~log(numingroup[i]) particles.
    //we set the limit at 2*log(numingroup[i]) to account for overhead in producing tree and calculating new potential
    iunbindsizeflag=(nEplus<2.0*log((double)nig));
    if (iunbindsizeflag==0) Potential(opt, nig, groupPart);
    else {
        for (auto k=0;k<nEplus;k++) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(r2,poti)
{
#pragma omp for schedule(dynamic) nowait
#endif
            for (auto j=0;j<nig;j++) {
                if (j!=nEplusid[k]) {
                    r2=0.;for (auto n=0;n<3;n++) r2+=pow(groupPart[nEplusid[k]].GetPosition(n)-groupPart[j].GetPosition(n),2.0);
                    r2+=eps2;
                    r2=1.0/sqrt(r2);
                    poti=groupPart[j].GetPotential()+opt.G*(groupPart[nEplusid[k]].GetMass()*groupPart[j].GetMass())*r2;
#ifdef NOMASS
                    poti*=mv2;
#endif
                    groupPart[j].SetPotential(poti);
                }
            }
#ifdef USEOPENMP
}
#endif
        }
    }
}

inline void RemoveGroup(Options &opt, Int_t &ning, Int_t *&pfof, Particle *&groupPart, int &iunbindflag)
{
    //if group too small remove entirely
    if (ning<opt.MinSize) {
        for (auto j=0;j<ning;j++) {
            pfof[groupPart[j].GetPID()]=0;
            groupPart[j].SetPID(-1);
        }
        ning=0;
        iunbindflag++;
    }
}

inline void AdjustPGListForUnbinding(int &unbindloops, Int_t &ning, Int_t *&pglist, Particle *&groupPart){
    if (unbindloops>0) for (auto j=0;j<ning;j++) pglist[j]=groupPart[j].GetPID();
}

inline void GetBoundFractionAndMaxE(Options &opt,
    Int_t &ning, Particle *groupPart, Coordinate &cmvel,
    Double_t &Efrac, Double_t &maxE, Int_t &nunbound,
    bool sortflag=true
    )
{
    Double_t mass, Ti, v2, totT=0;
    int nthreads;
    Efrac=nunbound=0;
    maxE=-MAXVALUE;
#ifdef USEOPENMP
    nthreads = max((int)(ning/(float)ompunbindnum),1);
    nthreads = min(nthreads,omp_get_max_threads());
    #pragma omp parallel for \
    default(shared) private(v2,Ti,mass) schedule(dynamic) \
    reduction(+:totT,Efrac,nunbound) reduction(max:maxE) num_threads(nthreads)
#endif
    for (auto j=0;j<ning;j++) {
        mass = groupPart[j].GetMass();
#ifdef NOMASS
        mass = opt.MassValue;
#endif
        v2=0.0;for (auto k=0;k<3;k++) v2+=pow(groupPart[j].GetVelocity(k)-cmvel[k],2.0);
        Ti=0.5*mass*v2;
#ifdef GASON
        Ti+=mass*groupPart[j].GetU();
#endif
        totT+=Ti;
        groupPart[j].SetDensity(opt.uinfo.Eratio*Ti+groupPart[j].GetPotential());
        if (groupPart[j].GetDensity()>maxE) maxE = groupPart[j].GetDensity();
        Efrac+=(Ti+groupPart[j].GetPotential()<0);
        if (groupPart[j].GetDensity()>0) nunbound++;
    }
    Efrac/=(Double_t)ning;
    //if object is not mostly unbound, then sort as will iteratively unbind
    if (nunbound<opt.uinfo.maxunboundfracforiterativeunbind*ning && sortflag) {
        sort(groupPart, groupPart + ning, [](const Particle &a, const Particle &b) {
        return a.GetDensity() < b.GetDensity();});
    }
}

//@}

///\name Remove unbound particles from a candidate group
//@{
/*!
    Interface for unbinding proceedure. Unbinding routine requires several arrays, such as numingroup, pglist,gPart,ids, etc
    This arrays may have been constructed prior to the unbinding call and so can be passed to the routine
    if this is called it uses Particle array then deletes it.
*/
int CheckUnboundGroups(Options opt, const Int_t nbodies, Particle *Part, Int_t &ngroup, Int_t *&pfof, Int_t *numingroup, Int_t **pglist, int ireorder, Int_t *groupflag)
{
    bool ningflag=false, pglistflag=false;
    int iflag;
    Int_t ng=ngroup;
    auto time1=MyGetTime();
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    if (opt.iverbose) cout<<ThisTask<<" Unbinding "<<ngroup<<" groups  ... "<<endl;

    //array creation check.
    if (numingroup==NULL) ningflag=true;
    if (pglist==NULL) pglistflag=true;
    //if array was not created outside, build now
    if (ningflag) numingroup=BuildNumInGroup(nbodies, ngroup, pfof);
    if (pglistflag) pglist=BuildPGList(nbodies, ngroup, numingroup, pfof);
    //unbind can either copy information into new arrays and not worry about order or reseting anything
    //OR save memory at the cost of a few more computations (specifically sorts).
    //Build particle array containing only particles in groups;
#ifdef SAVEMEM
    Int_t *storeval1=new Int_t[nbodies];
    Int_t *storeval2=new Int_t[nbodies];
    Int_t *noffset=new Int_t[ngroup+1];
    Double_t *storeden=new Double_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) {storeval1[i]=Part[i].GetPID();storeval2[i]=Part[i].GetID();storeden[i]=Part[i].GetDensity();Part[i].SetPID(i);Part[i].SetID(pfof[i]+nbodies*(pfof[i]==0));}
    // qsort(Part,nbodies,sizeof(Particle),IDCompare);
    std::sort(Part, Part + nbodies, IDCompareVec);
    noffset[0]=noffset[1]=0;
    for (Int_t i=2;i<=ngroup;i++) noffset[i]=noffset[i-1]+numingroup[i-1];
#else
    Particle **gPart=BuildPartList(ngroup, numingroup, pglist, Part);
    for (Int_t i=1;i<=ngroup;i++) {
        for (Int_t j=0;j<numingroup[i];j++) {
            gPart[i][j].SetID(j);
            gPart[i][j].SetPID(pglist[i][j]);
        }
    }
#endif

#ifdef SAVEMEM
    iflag=Unbind(opt, Part, ngroup, numingroup,noffset,pfof);
    //if keeping track of a flag, set flag to 0 if group no longer present
    if (groupflag!=NULL) {
        for (Int_t i=0;i<=ng;i++) if (numingroup[i]==0) groupflag[i]=0;
        //and reorder it if reordering group ids
        if (ireorder) ReorderGroupIDsAndArraybyValue(const Int_t numgroups, const Int_t newnumgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, Int_t *value, Int_t *gdata)
    }
    //reset order
    // qsort(Part,nbodies,sizeof(Particle),PIDCompare);
    std::sort(Part, Part + nbodies, PIDCompareVec);
    for (Int_t i=0;i<nbodies;i++) {Part[i].SetPID(storeval1[i]);Part[i].SetID(storeval2[i]);Part[i].SetDensity(storeden[i]);}
    //and alter pglist
    for (Int_t i=0;i<=ng;i++) numingroup[i]=0;
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) pglist[pfof[i]][numingroup[pfof[i]]++]=i;
    if (ireorder==1 && iflag&&ngroup>0) {
        if (groupflag!=NULL) {
            //allocate memory to store the group value which is a combination of the groupflag value and numingroup
            Double_t *groupvalue=new Double_t[ng+1];
            for (Int_t i=0;i<=ng;i++) {
                if (numingroup[i]==0) {groupvalue[i]=groupflag[i]=0;}
                else {
                    if (groupflag[i]==1) groupvalue[i]=numingroup[i];
                    else groupvalue[i]=1.0/(Double_t)numingroup[i];
                }
            }
            //and reorder it if reordering group ids
            if (ireorder) ReorderGroupIDsAndArraybyValue(ng,ngroup,numingroup,pfof,pglist,groupvalue,groupflag);
            delete[] groupvalue;
        }
        else ReorderGroupIDs(ng,ngroup,numingroup,pfof,pglist);
    }
    delete[] noffset;
#else

    //if groupflags are provided then explicitly reorder here if required, otherwise internal reordering within unbind.
    if (groupflag!=NULL) iflag = Unbind(opt, gPart, ngroup, numingroup,pfof,pglist,0);
    else iflag = Unbind(opt, gPart, ngroup, numingroup,pfof,pglist,ireorder);

    //if keeping track of a flag, set flag to 0 if group no longer present
    if (ireorder==1 && iflag&&ngroup>0) {
        if (groupflag!=NULL) {
            //allocate memory to store the group value which is a combination of the groupflag value and numingroup
            Double_t *groupvalue=new Double_t[ng+1];
            for (Int_t i=0;i<=ng;i++) {
                if (numingroup[i]==0) {groupvalue[i]=groupflag[i]=0;}
                else {
                    if (groupflag[i]==1) groupvalue[i]=numingroup[i];
                    else groupvalue[i]=1.0/(Double_t)numingroup[i];
                }
            }
            //and reorder it if reordering group ids
            if (ireorder) {
                ReorderGroupIDsAndArraybyValue(ng,ngroup,numingroup,pfof,pglist,groupvalue,groupflag);
            }
            delete[] groupvalue;
        }
        else {
            ReorderGroupIDs(ng,ngroup,numingroup,pfof,pglist);
        }
    }
    for (Int_t i=1;i<=ng;i++) delete[] gPart[i];delete[] gPart;
#endif
    if (pglistflag) {for (Int_t i=1;i<=ng;i++) delete[] pglist[i];delete[] pglist;}
    if (ningflag) delete[] numingroup;

    if (opt.iverbose) cout<<ThisTask<<" Done. Number of groups remaining "<<ngroup<<" in"<<MyElapsedTime(time1)<<endl;

    return iflag;
}

///Calculate potential of groups
inline void CalculatePotentials(Options &opt, Particle **gPart, Int_t &numgroups, Int_t *numingroup)
{
    int maxnthreads,nthreads=1;

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    if (!opt.uinfo.icalculatepotential) return;

    //for parallel environment store maximum number of threads
    nthreads=1;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=nthreads=omp_get_num_threads();
    }
#endif

    //for each group calculate potential
    //if group is small calculate potentials using PP
    //here openmp is over groups since each group is small
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (auto i=1;i<=numgroups;i++)
    {
        if (numingroup[i]<=POTPPCALCNUM) {
            if (numingroup[i]<0) continue;
            PotentialPP(opt, numingroup[i], gPart[i]);
        }
    }
#ifdef USEOPENMP
}
#endif
    //reset number of threads to maximum number
#ifdef USEOPENMP
#pragma omp master
    {
        omp_set_num_threads(maxnthreads);
    }
    nthreads=maxnthreads;
#endif
    for (auto i=1;i<=numgroups;i++)
    {
        if (numingroup[i]<0) continue;
        if (numingroup[i]>POTPPCALCNUM) {
            Potential(opt, numingroup[i], gPart[i]);
        }
    }
}

///Calculate potential of groups, assumes particle list is ordered by group
///and accessed by numingroup and noffset;
inline void CalculatePotentials(Options &opt, Particle *gPart, Int_t &numgroups, Int_t *numingroup, Int_t *noffset)
{
    int maxnthreads,nthreads=1;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    if (!opt.uinfo.icalculatepotential) return;

    //for parallel environment store maximum number of threads
    nthreads=1;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=nthreads=omp_get_num_threads();
    }
#endif

    //for each group calculate potential
    //if group is small calculate potentials using PP
    //here openmp is over groups since each group is small
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (auto i=1;i<=numgroups;i++)
    {
        if (numingroup[i]<=POTPPCALCNUM) {
            if (numingroup[i]<0) continue;
            PotentialPP(opt, numingroup[i], &gPart[noffset[i]]);
        }
    }
#ifdef USEOPENMP
}
#endif
    //reset number of threads to maximum number
#ifdef USEOPENMP
#pragma omp master
    {
        omp_set_num_threads(maxnthreads);
    }
    nthreads=maxnthreads;
#endif
    for (auto i=1;i<=numgroups;i++)
    {
        if (numingroup[i]<0) continue;
        if (numingroup[i]>POTPPCALCNUM) {
            Potential(opt, numingroup[i], &gPart[noffset[i]]);
        }
    }
}

///loop over groups and get velocity frame
inline void CalculateBindingReferenceFrame(Options &opt,
    Particle **gPart, Int_t &numgroups, Int_t *numingroup,
    Double_t *&gmass, Coordinate *&cmvel)
{
    Double_t potmin,menc;
    Int_t npot,ipotmin;
    Coordinate potpos;
    Int_t *storeval;

    //if using standard frame, then using CMVEL of the entire structure
    if (opt.uinfo.fracpotref==1.0) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
            for (auto i=1;i<=numgroups;i++)
            {
                for (auto k=0;k<3;k++) cmvel[i][k]=0;
                for (auto j=0;j<numingroup[i];j++) {
                    gmass[i]+=gPart[i][j].GetMass();
                    for (auto k=0;k<3;k++) cmvel[i][k]+=gPart[i][j].GetVelocity(k)*gPart[i][j].GetMass();
                }
                cmvel[i] *= 1.0/gmass[i];
            }
#ifdef USEOPENMP
}
#endif
    }
    else {
        if (opt.uinfo.cmvelreftype==CMVELREF) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(npot,menc,potpos,storeval)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
            for (auto i=1;i<=numgroups;i++)
            {
                if (numingroup[i]<0) continue;
                for (auto k=0;k<3;k++) potpos[k]=cmvel[i][k]=0;
                for (auto j=0;j<numingroup[i];j++) {
                    gmass[i]+=gPart[i][j].GetMass();
                    for (auto k=0;k<3;k++) potpos[k]+=gPart[i][j].GetPosition(k)*gPart[i][j].GetMass();
                }
                potpos*=(1.0/gmass[i]);
                //now sort by radius, first store original position in array
                storeval=new Int_t[numingroup[i]];
                for (auto j=0;j<numingroup[i];j++) {
                    for (auto k=0;k<3;k++) gPart[i][j].SetPosition(k,gPart[i][j].GetPosition(k)-potpos[k]);
                    storeval[i] = gPart[i][j].GetID();
                    gPart[i][j].SetID(j);
                }
                //sort by radius
                gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),RadCompare);
                //use central regions to define centre of mass velocity
                //determine how many particles to use
                npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
                npot=min(npot,numingroup[i]);
                cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
                for (auto j=0;j<npot;j++) {
                    for (auto k=0;k<3;k++) cmvel[i][k]+=gPart[i][j].GetVelocity(k)*gPart[i][j].GetMass();
                    menc+=gPart[i][j].GetMass();
                }
                for (auto j=0;j<3;j++) {cmvel[i][j]/=menc;}
                gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),IDCompare);
                for (auto j=0;j<numingroup[i];j++) {
                    gPart[i][j].SetID(storeval[j]);
                    for (auto k=0;k<3;k++) gPart[i][j].SetPosition(k,gPart[i][j].GetPosition(k)+potpos[k]);
                }
                delete[] storeval;
            }
#ifdef USEOPENMP
}
#endif
        }
        //if using potential then must identify minimum potential.
        //Note that  most computations involve sorts, so parallize over groups
        else if (opt.uinfo.cmvelreftype==POTREF) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(npot,menc,potmin,ipotmin,potpos,storeval)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
            for (auto i=1;i<=numgroups;i++) {
                if (numingroup[i]<0) continue;
                //determine how many particles to use
                npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
                npot=min(npot,numingroup[i]);

                storeval=new Int_t[numingroup[i]];
                for (auto j=0;j<numingroup[i];j++) {storeval[j]=gPart[i][j].GetID();gPart[i][j].SetID(j);}
                //determine position of minimum potential and by radius around this position
                potmin=gPart[i][0].GetPotential();ipotmin=0;
                for (auto j=1;j<numingroup[i];j++) if (gPart[i][j].GetPotential()<potmin) {potmin=gPart[i][j].GetPotential();ipotmin=j;}
                for (auto k=0;k<3;k++) potpos[k]=gPart[i][ipotmin].GetPosition(k);

                for (auto j=0;j<numingroup[i];j++) {
                    for (auto k=0;k<3;k++) gPart[i][j].SetPosition(k,gPart[i][j].GetPosition(k)-potpos[k]);
                }
                gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),RadCompare);
                //now determine kinetic frame
                cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
                for (auto j=0;j<npot;j++) {
                    for (auto k=0;k<3;k++) cmvel[i][k]+=gPart[i][j].GetVelocity(k)*gPart[i][j].GetMass();
                    menc+=gPart[i][j].GetMass();
                }
                for (auto j=0;j<3;j++) {cmvel[i][j]/=menc;}
                gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),IDCompare);
                for (auto j=0;j<numingroup[i];j++) gPart[i][j].SetID(storeval[j]);
                delete[] storeval;
            }
#ifdef USEOPENMP
}
#endif
        }
    }
}

///loop over groups and get velocity frame
inline void CalculateBindingReferenceFrame(Options &opt,
    Particle *gPart, Int_t &numgroups, Int_t *numingroup, Int_t *noffset,
    Double_t *&gmass, Coordinate *&cmvel)
{
    Double_t potmin,menc;
    Int_t npot,ipotmin;
    Coordinate potpos;
    Int_t *storeval;
    //if using standard frame, then using CMVEL of the entire structure
    if (opt.uinfo.fracpotref==1.0) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
            for (auto i=1;i<=numgroups;i++)
            {
                for (auto k=0;k<3;k++) cmvel[i][k]=0;
                for (auto j=0;j<numingroup[i];j++) {
                    gmass[i]+=gPart[noffset[i]+j].GetMass();
                    for (auto k=0;k<3;k++) cmvel[i][k]+=gPart[noffset[i]+j].GetVelocity(k)*gPart[noffset[i]+j].GetMass();
                }
                cmvel[i] *= 1.0/gmass[i];
            }
#ifdef USEOPENMP
}
#endif
    }
    else {
        if (opt.uinfo.cmvelreftype==CMVELREF) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(npot,menc,potpos,storeval)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
            for (auto i=1;i<=numgroups;i++)
            {
                if (numingroup[i]<0) continue;
                for (auto k=0;k<3;k++) potpos[k]=cmvel[i][k]=0;
                for (auto j=0;j<numingroup[i];j++) {
                    gmass[i]+=gPart[noffset[i]+j].GetMass();
                    for (auto k=0;k<3;k++) potpos[k]+=gPart[noffset[i]+j].GetPosition(k)*gPart[noffset[i]+j].GetMass();
                }
                potpos*=(1.0/gmass[i]);
                //now sort by radius, first store original position in array
                storeval=new Int_t[numingroup[i]];
                for (auto j=0;j<numingroup[i];j++) {
                    for (auto k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)-potpos[k]);
                    storeval[i] = gPart[noffset[i]+j].GetID();
                    gPart[noffset[i]+j].SetID(j);
                }
                //sort by radius
                gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),RadCompare);
                //use central regions to define centre of mass velocity
                //determine how many particles to use
                npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
                npot=min(npot,numingroup[i]);
                cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
                for (auto j=0;j<npot;j++) {
                    for (auto k=0;k<3;k++) cmvel[i][k]+=gPart[noffset[i]+j].GetVelocity(k)*gPart[noffset[i]+j].GetMass();
                    menc+=gPart[noffset[i]+j].GetMass();
                }
                for (auto j=0;j<3;j++) {cmvel[i][j]/=menc;}
                gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),IDCompare);
                for (auto j=0;j<numingroup[i];j++) {
                    gPart[noffset[i]+j].SetID(storeval[j]);
                    for (auto k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)+potpos[k]);
                }
                delete[] storeval;
            }
#ifdef USEOPENMP
}
#endif
        }
        //if using potential then must identify minimum potential.
        //Note that  most computations involve sorts, so parallize over groups
        else if (opt.uinfo.cmvelreftype==POTREF) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(npot,menc,potmin,ipotmin,potpos,storeval)
{
#pragma omp for schedule(dynamic,1) nowait
#endif
            for (auto i=1;i<=numgroups;i++) {
                if (numingroup[i]<0) continue;
                //determine how many particles to use
                npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
                npot=min(npot,numingroup[i]);

                storeval=new Int_t[numingroup[i]];
                for (auto j=0;j<numingroup[i];j++) {storeval[j]=gPart[noffset[i]+j].GetID();gPart[noffset[i]+j].SetID(j);}
                //determine position of minimum potential and by radius around this position
                potmin=gPart[noffset[i]+0].GetPotential();ipotmin=0;
                for (auto j=1;j<numingroup[i];j++) if (gPart[noffset[i]+j].GetPotential()<potmin) {potmin=gPart[noffset[i]+j].GetPotential();ipotmin=j;}
                for (auto k=0;k<3;k++) potpos[k]=gPart[noffset[i]+ipotmin].GetPosition(k);

                for (auto j=0;j<numingroup[i];j++) {
                    for (auto k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)-potpos[k]);
                }
                gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),RadCompare);
                //now determine kinetic frame
                cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
                for (auto j=0;j<npot;j++) {
                    for (auto k=0;k<3;k++) cmvel[i][k]+=gPart[noffset[i]+j].GetVelocity(k)*gPart[noffset[i]+j].GetMass();
                    menc+=gPart[noffset[i]+j].GetMass();
                }
                for (auto j=0;j<3;j++) {cmvel[i][j]/=menc;}
                gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),IDCompare);
                for (auto j=0;j<numingroup[i];j++) gPart[noffset[i]+j].SetID(storeval[j]);
                delete[] storeval;
            }
#ifdef USEOPENMP
}
#endif
        }
    }
}


/*!
    Unbinding algorithm that checks to see if a group is self-bound. For small groups the potential is calculated using a PP algorithm, for large groups a tree-potential using kd-tree and monopole is calculated. \n
    There are several ways a group can be defined as bound, it total energy is negative, its least bound particle has negative energy, or both. Also one can have different
    kinetic reference frames. By default, uses the CM velocity of the total system BUT could also use the velocity of a region centred on the minimum potential well. \n

    If a group is not self bound, the least bound particle is removed from the group. Then
    \arg if CM reference frame, the centre-of-mass velocity is recalculated and kinetic energies are recalculated (if it has changed enough). \n
    \arg the potential can be left unchanged or if one ignores the background unbound but linked particles, the energies are adjusted for the loss of the particle from the group.
    NOTE that for groups where the tree-potential was calculated, at the moment, the subtracted energy corresponds to the PP calculation which can lead to decrepancies. However, unless one is
    worried about the exact details of when an object is self-bound, this is not an issue. \n

    Finally, this routines assumes that the pglist passed to the routine is for a gPart array that was build in id order from pfof and a local particle array.
*/
int Unbind(Options &opt, Particle **gPart, Int_t &numgroups, Int_t *numingroup, Int_t *pfof, Int_t **pglist, int ireorder)
{
    //flag which is changed if any groups are altered as groups may need to be reordered.
    int iunbindflag=0;
    //flag used to determine what style of update to the potential is done for larger groups as
    //if the amount of particles removed is large enough for large groups, it is more efficient to
    //recalculate the entire potential using a Tree code than it is removing the contribution of each removed particle from
    //all other particles
    int maxnthreads,nthreads=1,n;
    Int_t i,j,k,ng=numgroups, oldnumingroup;
    int unbindloops;
    bool sortflag;
    Double_t maxE,v2,r2,poti,Ti, Efrac;
#ifdef NOMASS
    Double_t mv2=opt.MassValue*opt.MassValue;
#endif
    Int_t nEplus,maxunbindsize, nEfrac, nunbound;
    Int_t *nEplusid;
    int *Eplusflag;
    bool unbindcheck;

    Double_t *gmass;
    Coordinate *cmvel;

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    //note that it is possible that called as a library, velociraptor
    //does not need to calculate potentials itself
    //in that case do not calculate potentials but instead
    //copy relevant information
    cmvel   =new Coordinate[numgroups+1];
    gmass   =new Double_t[numgroups+1];
    for (i=1;i<=numgroups;i++) {
        cmvel[i]=Coordinate(0.);
        gmass[i]=0.;
        if (opt.uinfo.icalculatepotential) {
            for (j=0;j<numingroup[i];j++) gPart[i][j].SetPotential(0);
        }
        #ifdef SWIFTINTERFACE
        else {
            for (j=0;j<numingroup[i];j++) gPart[i][j].SetPotential(gPart[i][j].GetGravityPotential());
        }
        #endif
    }

    //if calculate potential
    CalculatePotentials(opt, gPart, numgroups, numingroup);

    //Now set the kinetic reference frame
    CalculateBindingReferenceFrame(opt, gPart, numgroups, numingroup, gmass, cmvel);

    //now go through groups and begin unbinding by finding least bound particle
    //again for small groups multithread over groups
    //larger groups thread over particles in a group
    //for large groups, paralleize over particle, for small groups parallelize over groups
    //here energy data is stored in density
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,n,maxE,maxunbindsize,nEplus,nEplusid,Eplusflag,v2,Ti,unbindcheck,Efrac,nEfrac,nunbound,r2,poti,unbindloops,sortflag,oldnumingroup)
{
    #pragma omp for schedule(dynamic) nowait reduction(+:iunbindflag)
#endif
    for (i=1;i<=numgroups;i++) if (numingroup[i]<ompunbindnum && numingroup[i]>0)
    {
        unbindloops=0;
        oldnumingroup = numingroup[i];
        GetBoundFractionAndMaxE(opt, numingroup[i], gPart[i], cmvel[i], Efrac, maxE,nunbound);
        if (nunbound>=opt.uinfo.maxunboundfracforiterativeunbind*numingroup[i]) {
            for (j=0;j<numingroup[i];j++) pfof[pglist[i][j]]=0;
            numingroup[i]=0;
            iunbindflag++;
        }
        else {
            //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
            maxunbindsize=(Int_t)(opt.uinfo.maxunbindfrac*nunbound+1);
            Eplusflag=new int[numingroup[i]];
            nEplusid=new Int_t[numingroup[i]];
            unbindcheck = CheckGroupForBoundness(opt,Efrac,maxE,numingroup[i]);
            FillUnboundArrays(opt, maxunbindsize, numingroup[i], gPart[i], Efrac, nEplusid, Eplusflag, nEplus, unbindcheck);
            while(unbindcheck)
            {
                iunbindflag++;
                unbindloops++;
                UpdateCMForUnboundParticles(opt, gmass[i], cmvel[i],
                    numingroup[i], gPart[i], nEplus, nEplusid, Eplusflag);
                UpdatePotentialForUnboundParticlesPP(opt, numingroup[i], gPart[i],
                    nEplus, nEplusid, Eplusflag);
                //remove particles with positive energy
                RemoveUnboundParticles(i, pfof, numingroup[i], pglist[i], gPart[i], nEplus, nEplusid, Eplusflag);
                //if number of particles remove with positive energy is near to the number allowed to be removed
                //must recalculate kinetic energies and check if maxE>0
                //otherwise, end unbinding.
                if (nEplus<opt.uinfo.maxallowedunboundfrac*numingroup[i]) {
                    unbindcheck=false;
                }
                else {
                    sortflag=false;
                    if ((oldnumingroup-numingroup[i])>opt.uinfo.maxallowedunboundfrac*oldnumingroup) {
                        oldnumingroup=numingroup[i];
                        sortflag=true;
                    }
                    //recalculate kinetic energies since cmvel has changed
                    GetBoundFractionAndMaxE(opt, numingroup[i], gPart[i], cmvel[i], Efrac, maxE,nunbound, sortflag);
                    maxunbindsize=(Int_t)(opt.uinfo.maxunbindfrac*nunbound+1);
                    unbindcheck = CheckGroupForBoundness(opt,Efrac,maxE,numingroup[i]);
                    FillUnboundArrays(opt, maxunbindsize, numingroup[i], gPart[i], Efrac, nEplusid, Eplusflag, nEplus, unbindcheck);
                }
            }
            //if group too small remove entirely
            AdjustPGListForUnbinding(unbindloops,numingroup[i],pglist[i],gPart[i]);
            RemoveGroup(opt, numingroup[i], pfof, gPart[i], iunbindflag);
            delete[] nEplusid;
            delete[] Eplusflag;
        }
    }
#ifdef USEOPENMP
}
#endif
    for (i=1;i<=numgroups;i++) if (numingroup[i]>=ompunbindnum)
    {
        unbindloops=0;
        oldnumingroup = numingroup[i];
        GetBoundFractionAndMaxE(opt, numingroup[i], gPart[i], cmvel[i], Efrac, maxE, nunbound);
        //if amount unbound is very large, just remove group entirely
        if (nunbound>=opt.uinfo.maxunboundfracforiterativeunbind*numingroup[i]) {
            for (j=0;j<numingroup[i];j++) pfof[pglist[i][j]]=0;
            numingroup[i]=0;
            iunbindflag++;
        }
        else {
            //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
            maxunbindsize=(Int_t)(opt.uinfo.maxunbindfrac*nunbound+1);
            nEplusid=new Int_t[numingroup[i]];
            Eplusflag=new int[numingroup[i]];
            //check if bound;
            unbindcheck = CheckGroupForBoundness(opt,Efrac,maxE,numingroup[i]);
            FillUnboundArrays(opt, maxunbindsize, numingroup[i], gPart[i], Efrac, nEplusid, Eplusflag, nEplus, unbindcheck);
            while(unbindcheck)
            {
                iunbindflag+=1;
                unbindloops++;
                UpdateCMForUnboundParticles(opt, gmass[i], cmvel[i],
                    numingroup[i], gPart[i], nEplus, nEplusid, Eplusflag);
                UpdatePotentialForUnboundParticles(opt, numingroup[i], gPart[i],
                    nEplus, nEplusid, Eplusflag);
                RemoveUnboundParticles(i, pfof, numingroup[i], pglist[i], gPart[i], nEplus, nEplusid, Eplusflag);
                //if number of particles remove with positive energy is near to the number allowed to be removed
                //must recalculate kinetic energies and check if maxE>0
                //otherwise, end unbinding.
                if (nEplus<opt.uinfo.maxallowedunboundfrac*numingroup[i]) {
                    unbindcheck=false;
                    continue;
                }
                else{
                    sortflag=false;
                    if ((oldnumingroup-numingroup[i])>opt.uinfo.maxallowedunboundfrac*oldnumingroup) {
                        oldnumingroup=numingroup[i];
                        sortflag=true;
                    }
                    //recalculate kinetic energies since cmvel has changed
                    GetBoundFractionAndMaxE(opt, numingroup[i], gPart[i], cmvel[i], Efrac, maxE, nunbound,sortflag);
                    //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
                    maxunbindsize=(Int_t)(opt.uinfo.maxunbindfrac*nunbound+1);
                    unbindcheck = CheckGroupForBoundness(opt,Efrac,maxE,numingroup[i]);
                    FillUnboundArrays(opt, maxunbindsize, numingroup[i], gPart[i], Efrac, nEplusid, Eplusflag, nEplus, unbindcheck);
                }
            }
            AdjustPGListForUnbinding(unbindloops,numingroup[i],pglist[i],gPart[i]);
            RemoveGroup(opt, numingroup[i], pfof, gPart[i], iunbindflag);
            delete[] nEplusid;
            delete[] Eplusflag;
        }
    }

    for (i=1;i<=numgroups;i++) if (numingroup[i]==0) ng--;
    if (ireorder==1 && iunbindflag&&ng>0) ReorderGroupIDs(numgroups,ng,numingroup,pfof,pglist);
    delete[] cmvel;
    delete[] gmass;
    numgroups=ng;
    //return if any unbinding done indicating groups have been reordered
    if (iunbindflag) return 1;
    else return 0;
}

/// Calculates the gravitational potential using a kd-tree and monopole expansion
///\todo need ewald correction for periodic systems and also use more than monopole.
void Potential(Options &opt, Int_t nbodies, Particle *Part, Double_t *potV)
{
    Potential(opt, nbodies, Part);
    for (auto i=0;i<nbodies;i++) potV[i]=Part[i].GetPotential();
}

void Potential(Options &opt, Int_t nbodies, Particle *Part)
{
    Int_t oldnbodies;
    KDTree *tree;
    Particle *part;
    bool runomp = false;
    int nsearch;
    double mr;
    int bsize = opt.uinfo.BucketSize;

    ///\todo need to get nomass stuff working

    //if approximate potential calculated, subsample partile distribution
    oldnbodies = nbodies;
    ParticleSubSample(opt, oldnbodies, Part, nbodies, part, mr);
    if (part != Part) bsize = ceil(bsize*opt.uinfo.approxpotnumfrac);

    //build tree and calculate a tree based potential
    tree = new KDTree(part, nbodies, bsize, tree->TPHYS, tree->KEPAN,
        100, 0, 0, 0, NULL, NULL, runomp);
    if (part != Part) tree->OverWriteInputOrder();
    PotentialTree(opt, nbodies, part, tree);
    //and assign potentials back if running approximate potential calculation
    //i.e., particle pointer does not point to original particle pointer
    if (part != Part) {
        nsearch = min(4,(int)ceil(mr+1));
        PotentialInterpolate(opt, oldnbodies, Part, part, tree, mr, nsearch);
    }
    delete tree;

    //if iterating then run again but new shuffled indices
    //to better sample low potential regions
    // if (opt.uinfo.approxpotiterate > 0 && part != Part && nbodies<0.5*oldnbodies) {
    //     int numloops = 0;
    //     vector<indexpot> potindex(oldnbodies);
    //     double newpotsum, oldpotsum, deltapot, oldminpot, newminpot, deltaminpot;
    //     double oldsubpotsum, newsubpotsum;
    //     // double deltapot2;
    //     oldpotsum = 0;
    //     oldminpot = newminpot = 0;
    //     for (auto i=0;i<oldnbodies;i++) {
    //         potindex[i].potval = Part[i].GetPotential();
    //         potindex[i].index = i;
    //         oldpotsum += potindex[i].potval;
    //         if (potindex[i].potval<oldminpot) oldminpot = potindex[i].potval;
    //     }
    //     oldsubpotsum = 0; for (auto i=0;i<nbodies;i++) oldsubpotsum += part[i].GetPotential();
    //
    //     cout<<oldnbodies<<" "<<nbodies<<" beginning iteration "<<oldpotsum<<endl;
    //     do {
    //         //sort indices by pot value
    //         sort(potindex.begin(), potindex.end(), [](indexpot &a, indexpot &b){
    //         return a.potval < b.potval;
    //         });
    //         //take random sample of the region around the min potential
    //         random_shuffle(potindex.begin(), potindex.begin()+oldnbodies/2);
    //         for (auto i=0;i<nbodies/2;i++) {
    //             Int_t index = potindex[i].index;
    //             part[i] = Part[index];
    //             part[i].SetMass(part[i].GetMass()*mr);
    //         }
    //         //and outer region
    //         random_shuffle(potindex.begin()+oldnbodies/2, potindex.end());
    //         for (auto i=nbodies/2;i<nbodies;i++) {
    //             Int_t index = potindex[i+oldnbodies/2].index;
    //             part[i] = Part[index];
    //             part[i].SetMass(part[i].GetMass()*mr);
    //         }
    //
    //         //build tree and calculate a tree based potential
    //         tree = new KDTree(part, nbodies, bsize, tree->TPHYS, tree->KEPAN,
    //             100, 0, 0, 0, NULL, NULL, runomp);
    //         tree->OverWriteInputOrder();
    //         PotentialTree(opt, nbodies, part, tree);
    //         PotentialInterpolate(opt, oldnbodies, Part, part, tree, mr, nsearch);
    //         delete tree;
    //
    //         // for (auto i=0;i<nbodies/2;i++) {
    //         //     Int_t index = potindex[i].index;
    //         //     part[i] = Part[index];
    //         //     part[i].SetMass(part[i].GetMass());
    //         // }
    //         // //and outer region
    //         // random_shuffle(potindex.begin()+nbodies/2, potindex.end());
    //         // for (auto i=nbodies/2;i<nbodies;i++) {
    //         //     Int_t index = potindex[i].index;
    //         //     part[i] = Part[index];
    //         //     part[i].SetMass(part[i].GetMass()*mr);
    //         // }
    //         //
    //         // //build tree and calculate a tree based potential
    //         // tree = new KDTree(part, nbodies, bsize, tree->TPHYS, tree->KEPAN,
    //         //     100, 0, 0, 0, NULL, NULL, runomp);
    //         // tree->OverWriteInputOrder();
    //         // PotentialTree(opt, nbodies, part, tree);
    //         // Particle *pblah = &Part[nbodies/2];
    //         // PotentialInterpolate(opt, oldnbodies-nbodies/2, pblah, part, tree, mr, nsearch);
    //         // pblah = NULL;
    //         // delete tree;
    //         // for (auto i=0;i<nbodies/2;i++) {
    //         //     Int_t index = potindex[i].index;
    //         //     Part[index].SetPotential(potindex[i].potval);
    //         // }
    //         //
    //         // newsubpotsum = 0; for (auto i=0;i<nbodies;i++) newsubpotsum += part[i].GetPotential();
    //         // cout<<oldnbodies<<" "<<nbodies<<" sub pot "<<oldsubpotsum<<" "<<newsubpotsum<<endl;
    //
    //         newpotsum = 0;
    //         newminpot = 0;
    //         for (auto i=0;i<oldnbodies;i++) {
    //             potindex[i].potval = Part[i].GetPotential();
    //             potindex[i].index = i;
    //             newpotsum += potindex[i].potval;
    //             if (potindex[i].potval<newminpot) newminpot = potindex[i].potval;
    //         }
    //         deltapot = (oldpotsum-newpotsum)/newpotsum;
    //         deltapot = sqrt(deltapot*deltapot);
    //         deltaminpot = (oldminpot-newminpot)/newminpot;
    //         deltaminpot = sqrt(deltaminpot*deltaminpot);
    //         numloops++;
    //         // deltapot2 = (potfull-newpotsum)/potfull;
    //         // deltapot2 = sqrt(deltapot2*deltapot2);
    //         //
    //         cout<<oldnbodies<<" "<<nbodies<<" Iterating "<<oldpotsum<<" "<<newpotsum<<" "<<deltapot<<" "<<numloops<<" "<<deltaminpot<<endl;
    //         oldpotsum = newpotsum;
    //         oldminpot = newminpot;
    //     } while((deltapot>opt.uinfo.approxpotrelerr || deltaminpot>opt.uinfo.approxpotrelerr) && numloops<opt.uinfo.approxpotiterate);
    // }

    //free up memory
    if (part != Part) delete[] part;
    else part = NULL;
}

void ParticleSubSample(Options &opt, const Int_t nbodies, Particle *&Part,
    Int_t &newnbodies, Particle *&newpart, double &mr)
{
    //if approximate potential calculated, subsample partile distribution
    newnbodies = nbodies;
    newpart = Part;
    mr = 1.0;
    if (opt.uinfo.iapproxpot) {
        if (newnbodies > opt.uinfo.approxpotminnum) {
            newnbodies = max((Int_t)(opt.uinfo.approxpotnumfrac*nbodies), (Int_t)opt.uinfo.approxpotminnum);
        }
        if (newnbodies < 0.5*nbodies) {
            if (opt.uinfo.approxpotmethod == POTAPPROXMETHODTREE) {
                //build tree that contains leaf nodes containing the desired
                //number of particles per leaf node
                Int_t bsize = ceil(nbodies/(float)newnbodies);
                KDTree *tree;
                tree = new KDTree(Part, nbodies, bsize, tree->TPHYS,tree->KEPAN,100);
                //first get all local leaf nodes;
                newnbodies = tree->GetNumLeafNodes();
                newpart = new Particle[newnbodies];
                mr = (double)nbodies/(double)newnbodies;
                Node *node;
                vector<leaf_node_info> leafnodes(newnbodies);
                Int_t inode=0, ipart=0;
                while (ipart<nbodies) {
                    node=tree->FindLeafNode(ipart);
                    leafnodes[inode].id = inode;
                    leafnodes[inode].istart = node->GetStart();
                    leafnodes[inode].iend = node->GetEnd();
                    leafnodes[inode].numtot = node->GetCount();
                    ipart+=leafnodes[inode].numtot;
                    inode++;
                }
                node=NULL;

#ifdef USEOPENMP
#pragma omp parallel for default(shared) schedule(dynamic)
#endif
                for (auto i=0;i<newnbodies;i++)
                {
                    double mass = 0;
                    leafnodes[i].cm[0]=leafnodes[i].cm[1]=leafnodes[i].cm[2]=0;
                    for (auto j=leafnodes[i].istart;j<leafnodes[i].iend;j++)
                    {
                        mass += Part[j].GetMass();
                        for (auto k=0;k<3;k++) leafnodes[i].cm[k] += Part[j].GetPosition(k)*Part[j].GetMass();
                    }
                    for (auto k=0;k<3;k++) newpart[i].SetPosition(k, leafnodes[i].cm[k]/mass);
                    newpart[i].SetMass(mass);
                }
                leafnodes.clear();
                delete tree;
            }
            else if (opt.uinfo.approxpotmethod == POTAPPROXMETHODRAND) {
                //randomly sample particle distribution
                mr = (double)nbodies/(double)newnbodies;
                newpart = new Particle[newnbodies];
                vector<Int_t> indices(nbodies);
                std::random_device rd;
                std::mt19937 g(rd());
                for (auto i=0;i<nbodies;i++) indices[i] = i;
                std::shuffle(indices.begin(), indices.end(), g);
                for (auto i=0;i<newnbodies;i++) {
                    Int_t index = indices[i];
                    newpart[i] = Part[index];
                    newpart[i].SetMass(newpart[i].GetMass()*mr);
                }
                indices.clear();
            }
        }
        else newnbodies = nbodies;
    }
}

void PotentialTree(Options &opt, Int_t nbodies, Particle *&Part, KDTree* &tree)
{
    Int_t ntreecell, nleafcell;
    Double_t r2, eps2=opt.uinfo.eps*opt.uinfo.eps;
#ifdef NOMASS
    Double_t mv2=opt.MassValue*opt.MassValue;
#endif
    int bsize = opt.uinfo.BucketSize;
    int maxnthreads, nthreads = 1;
    //for tree code potential calculation
    Int_t ncell;
    Int_t *start,*end;
    Double_t *cmtot, *cBmax,*cR2max;
    vector<Double_t>  r2val;
    Coordinate *cellcm;
    Node *root;
    Node **nodelist;
    Node* npomp;
    vector<Int_t> marktreecell, markleafcell;
    bool runomp = false;
#ifdef USEOPENMP
    runomp = (nbodies > POTOMPCALCNUM);
    #pragma omp parallel
    {
        if (omp_get_thread_num() == 0) maxnthreads = nthreads = omp_get_num_threads();
    }
#endif

    ncell=tree->GetNumNodes();
    root=tree->GetRoot();
    //to store particles in a given node
    start=new Int_t[ncell];
    end=new Int_t[ncell];
    //distance calculations used to determine when one uses cell or when one uses particles
    cmtot=new Double_t[ncell];
    cBmax=new Double_t[ncell];
    cR2max=new Double_t[ncell];
    cellcm=new Coordinate[ncell];
    //to store note list
    nodelist=new Node*[ncell];

    //from root node calculate cm for each node
    //start at root node and recursively move through list
    ncell=0;
    GetNodeList(root,ncell,nodelist,bsize);
    ncell++;

    //determine cm for all cells and openings
#ifdef USEOPENMP
#pragma omp parallel default(shared) if (runomp)
{
    #pragma omp for schedule(static)
#endif
    for (auto j=0;j<ncell;j++) {
        start[j]=(nodelist[j])->GetStart();
        end[j]=(nodelist[j])->GetEnd();
        cellcm[j][0]=cellcm[j][1]=cellcm[j][2]=0.;
        cmtot[j]=0;
        for (auto k=start[j];k<end[j];k++) {
            for (auto n=0;n<3;n++) cellcm[j][n]+=Part[k].GetPosition(n)*Part[k].GetMass();
            cmtot[j]+=Part[k].GetMass();
        }
        for (auto n=0;n<3;n++) cellcm[j][n]/=cmtot[j];
        Double_t xdiff,xdiff1;
        xdiff=(cellcm[j]-Coordinate(Part[start[j]].GetPosition())).Length();
        for (auto k=start[j]+1;k<end[j];k++) {
            xdiff1=(cellcm[j]-Coordinate(Part[k].GetPosition())).Length();
            if (xdiff<xdiff1) xdiff=xdiff1;
        }
        cBmax[j]=xdiff;
        cR2max[j]=4.0/3.0*xdiff*xdiff/(opt.uinfo.TreeThetaOpen*opt.uinfo.TreeThetaOpen);
    }
#ifdef USEOPENMP
}
#endif

    //then for each cell find all other cells that contain Particles within a cells gRmax and mark those
    //and mark all cells for which one does not have to unfold
    //for marked cells calculate pp, for every other cell just use the CM of the cell to calculate the potential.
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(marktreecell, markleafcell, r2val, npomp, ntreecell, nleafcell, r2) \
if (runomp)
{
#endif
    r2val.resize(ncell);
    marktreecell.resize(ncell);
    markleafcell.resize(ncell);
#ifdef USEOPENMP
    #pragma omp for schedule(static)
#endif
    for (auto j=0;j<nbodies;j++) {
        npomp=tree->GetRoot();
        Part[j].SetPotential(0.);
        ntreecell=nleafcell=0;
        Coordinate xpos(Part[j].GetPosition());
        MarkCell(npomp, marktreecell.data(), markleafcell.data(), 
            ntreecell, nleafcell, r2val.data(), 
            bsize, cR2max, cellcm, cmtot, xpos, eps2);
        for (auto k=0;k<ntreecell;k++) {
          Part[j].SetPotential(Part[j].GetPotential()-Part[j].GetMass()*r2val[k]);
        }
        for (auto k=0;k<nleafcell;k++) {
            for (auto l=start[markleafcell[k]];l<end[markleafcell[k]];l++) {
                if (j!=l) {
                    r2=0.;
                    for (auto n=0;n<3;n++) r2+=pow(Part[j].GetPosition(n)-Part[l].GetPosition(n),(Double_t)2.0);
                    r2+=eps2;
                    r2=1.0/sqrt(r2);
                    Part[j].SetPotential(Part[j].GetPotential()-(Part[j].GetMass()*Part[l].GetMass())*r2);
                }
            }
        }
        Part[j].SetPotential(Part[j].GetPotential()*opt.G);
#ifdef NOMASS
        Part[j].SetPotential(Part[j].GetPotential()*mv2);
#endif
    }
#ifdef USEOPENMP
}
#endif

    delete[] start;
    delete[] end;
    delete[] cmtot;
    delete[] cBmax;
    delete[] cR2max;
    delete[] cellcm;
    delete[] nodelist;
}

void PotentialInterpolate(Options &opt, const Int_t nbodies, Particle *&Part, Particle *&interpolatepart, KDTree *&tree, double massratio, int nsearch)
{
    bool runomp = false;
    vector<Int_t> nn;
    vector<Double_t> dist2;
    Double_t pot, wsum, w;
    runomp = (nbodies > POTOMPCALCNUM);
#ifdef USEOPENMP
#pragma omp parallel default(shared) private(nn, dist2, wsum, pot) \
if (runomp)
{
#endif
    nn.resize(nsearch);
    dist2.resize(nsearch);
#ifdef USEOPENMP
#pragma omp for schedule(static)
#endif
    for (auto i=0; i<nbodies; i++)
    {
        tree->FindNearestPos(Part[i].GetPosition(), nn.data(), dist2.data(), nsearch);
        pot = 0;
        wsum = 0;
        if (dist2[0] == 0) {
            pot = interpolatepart[nn[0]].GetPotential();
        }
        else {
            for (auto j=0;j<nsearch;j++) {
                w = 1.0/sqrt(dist2[j]);
                pot += interpolatepart[nn[j]].GetPotential()*w;
                wsum += w;
            }
            pot /= wsum;
        }
        Part[i].SetPotential(pot/massratio);
    }
    nn.clear();
    dist2.clear();
#ifdef USEOPENMP
}
#endif
}


void PotentialPP(Options &opt, Int_t nbodies, Particle *Part)
{
    Double_t r2, pot, poti, eps2=opt.uinfo.eps*opt.uinfo.eps;
#ifdef NOMASS
    Double_t mv2=opt.MassValue*opt.MassValue;
#endif
    for (auto j=0;j<nbodies;j++) Part[j].SetPotential(0.);
    for (auto j=0;j<nbodies;j++) {
        for (auto k=j+1;k<nbodies;k++) {
            r2=0.;for (auto n=0;n<3;n++) r2+=pow(Part[j].GetPosition(n)-Part[k].GetPosition(n),2.0);
            r2+=eps2;
            r2=1.0/sqrt(r2);
            pot=-opt.G*(Part[j].GetMass()*Part[k].GetMass())*r2;
            poti=Part[j].GetPotential()+pot;Part[j].SetPotential(poti);
            poti=Part[k].GetPotential()+pot;Part[k].SetPotential(poti);
        }
    }
    #ifdef NOMASS
    for (auto j=0;j<nbodies;j++) Part[j].SetPotential(Part[j].GetPotential()*mv2);
    #endif
}
