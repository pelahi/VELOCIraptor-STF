/*! \file binding.cxx
 *  \brief this file contains routines to calculate the binding energy
 *  Calculates gravitational energy using tree code, adds thermal self-energy if necessary for gas particles
 *  sorts, the particle list according to the energy and adds on a simple triaxial NFW potential based on
 *  the halo data describing the surrounding dm halo
 *  \todo need to add routines for ionizing radiation from stars, and possible other non-standard dm effects
 */

#include "baryoniccontent.h"

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

///calculate potential
void Potential(Options &opt, Int_t nbodies, Particle *Part)
{
    int maxnthreads,nthreads,l,n;
    Int_t i,j,k,ntreecell,nleafcell;
    Double_t v2,r2,Ti,eps2=opt.uinfo.eps*opt.uinfo.eps;
    //for tree code potential calculation
    Int_t ncell;
    Int_t *start,*end;
    Double_t *cmtot,*cBmax,*cR2max, **r2val;
    Coordinate *cellcm;
    Node *root;
    Node **nodelist, **npomp;
    Int_t **marktreecell,**markleafcell;
    //Double_t **nnr2;
    KDTree *tree;

    //for parallel environment store maximum number of threads
    nthreads=1;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=nthreads=omp_get_num_threads();
    }
#endif
    //otherwise use tree tree gravity calculation
    //here openmp is per group since each group is large
    //to make this memory efficient really need just KDTree that uses Coordinates
    tree=new KDTree(Part,nbodies,opt.uinfo.BucketSize,tree->TPHYS);
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

    //search tree
    marktreecell=new Int_t*[nthreads];
    markleafcell=new Int_t*[nthreads];
    r2val=new Double_t*[nthreads];
    npomp=new Node*[nthreads];
    for (j=0;j<nthreads;j++) {marktreecell[j]=new Int_t[ncell];markleafcell[j]=new Int_t[ncell];}
    for (j=0;j<nthreads;j++) {r2val[j]=new Double_t[ncell];}
    //from root node calculate cm for each node
    //start at root node and recursively move through list
    ncell=0;
    GetNodeList(root,ncell,nodelist,opt.uinfo.BucketSize);
    ncell++;
    //determine cm for all cells and openings
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,k,n)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (j=0;j<ncell;j++) {
        start[j]=(nodelist[j])->GetStart();
        end[j]=(nodelist[j])->GetEnd();
        cellcm[j][0]=cellcm[j][1]=cellcm[j][2]=0.;
        cmtot[j]=0;
        for (k=start[j];k<end[j];k++) {
            for (n=0;n<3;n++) cellcm[j][n]+=Part[k].GetPosition(n)*Part[k].GetMass();
            cmtot[j]+=Part[k].GetMass();
        }
        for (n=0;n<3;n++) cellcm[j][n]/=cmtot[j];
        Double_t xdiff,xdiff1;
        xdiff=(cellcm[j]-Coordinate(Part[start[j]].GetPosition())).Length();
        for (k=start[j]+1;k<end[j];k++) {
            xdiff1=(cellcm[j]-Coordinate(Part[k].GetPosition())).Length();
            if (xdiff<xdiff1) xdiff=xdiff1;
        }
        cBmax[j]=xdiff;
        cR2max[j]=4.0/3.0*xdiff*xdiff/(opt.uinfo.TreeThetaOpen*opt.uinfo.TreeThetaOpen);
    }
#ifdef USEOPENMP
}
#endif
    //then for each cell find all other cells that contain particles within a cells gRmax and mark those
    //and mark all cells for which one does not have to unfold
    //for marked cells calculate pp, for every other cell just use the CM of the cell to calculate the potential.
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,k,l,n,ntreecell,nleafcell,r2)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (j=0;j<nbodies;j++) {
        int tid;
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        npomp[tid]=tree->GetRoot();
        Part[j].SetPotential(0.);
        ntreecell=nleafcell=0;
        Coordinate xpos(Part[j].GetPosition());
        MarkCell(npomp[tid],marktreecell[tid], markleafcell[tid],ntreecell,nleafcell,r2val[tid],opt.uinfo.BucketSize, cR2max, cellcm, cmtot, xpos, eps2);
        for (k=0;k<ntreecell;k++) {
          Part[j].SetPotential(Part[j].GetPotential()-Part[j].GetMass()*r2val[tid][k]);
        }
        for (k=0;k<nleafcell;k++) {
            for (l=start[markleafcell[tid][k]];l<end[markleafcell[tid][k]];l++) {
                if (j!=l) {
                    r2=0.;for (n=0;n<3;n++) r2+=pow(Part[j].GetPosition(n)-Part[l].GetPosition(n),(Double_t)2.0);
                    r2+=eps2;
                    r2=1.0/sqrt(r2);
                    Part[j].SetPotential(Part[j].GetPotential()-(Part[j].GetMass()*Part[l].GetMass())*r2);
                }
            }
        }
        Part[j].SetPotential(Part[j].GetPotential()*opt.G);
    }
#ifdef USEOPENMP
}
#endif
    delete tree;
}

//@}

///Gas energy routines, assumes that particle class has temparture
/*
inline Double_t GetInternalEnergy(Double_t Boltzmann,Particle &p){
  return 1.5*(exp(Boltzmann+log(p.GetTemp()));
}
*/

///gets binding energy assuming particles are in reference frame and stores the ratio of kinetic to potential
///in GetPotential()
void GetBindingEnergy(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp) {
    cout<<"Get Energy"<<endl;
    Particle *Pval;
    Int_t i,j,k;
    Double_t eps2=opt.uinfo.eps*opt.uinfo.eps;
    Double_t r2,v2,Ti,poti,pot;
    Double_t Tval,Potval,Efracval,Eval,intE;
    Double_t *Tvaltyped,*Potvaltyped,*Efracvaltyped;
    Double_t mv2=opt.MassValue*opt.MassValue;
    Int_t noffset;
    Double_t t1=MyGetTime();
    Int_t ngdone=0;
    int maxnthreads=1,nthreads=1,tid;
    Int_t omparraysize;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=nthreads=omp_get_num_threads();
    }
#endif
    omparraysize=nthreads*NPARTTYPES;
    Tvaltyped=new Double_t[omparraysize];
    Potvaltyped=new Double_t[omparraysize];
    Efracvaltyped=new Double_t[omparraysize];

#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,r2,v2,poti,Ti,pot,Eval,intE,noffset)
{
    #pragma omp for schedule(dynamic,1) reduction(+:ngdone) nowait
#endif
    for (i=0;i<ngroup;i++) if (hp[i].AllNumberofParticles<ompunbindnum && hp[i].AllNumberofParticles>0) {
        noffset=hp[i].noffset;
        if (opt.ipotcalc) {
            pdata[i].T=pdata[i].Pot=pdata[i].Efrac=0.;
            for (j=0;j<NPARTTYPES;j++) pdata[i].Ttyped[j]=pdata[i].Efractyped[j]=0.;
        }
        else {
            pdata[i].T=pdata[i].Efrac=0.;
            for (j=0;j<NPARTTYPES;j++) pdata[i].Ttyped[j]=pdata[i].Efractyped[j]=0.;
        }
        if (opt.ipotcalc) {
        for (j=0;j<hp[i].AllNumberofParticles;j++) {
            for (k=j+1;k<hp[i].AllNumberofParticles;k++) {
                r2=0.;for (int n=0;n<3;n++) r2+=pow(Part[j+noffset].GetPosition(n)-Part[k+noffset].GetPosition(n),2.0);
                r2+=eps2;
                r2=1.0/sqrt(r2);
#ifdef NOMASS
                pot=-opt.G*(Part[j+noffset].GetMass()*Part[k+noffset].GetMass())*r2*mv2;
#else
                pot=-opt.G*(Part[j+noffset].GetMass()*Part[k+noffset].GetMass())*r2;
#endif
                poti=Part[j+noffset].GetPotential()+pot;Part[j+noffset].SetPotential(poti);
                poti=Part[k+noffset].GetPotential()+pot;Part[k+noffset].SetPotential(poti);
            }
            pdata[i].Pot+=Part[j+noffset].GetPotential();
            pdata[i].Pottyped[Part[j+noffset].GetType()]+=Part[j+noffset].GetPotential();
        }
        }
        for (j=0;j<hp[i].AllNumberofParticles;j++) {
            v2=0.;for (int n=0;n<3;n++) v2+=pow(Part[j+noffset].GetVelocity(n),2.0);
            intE=0;
            //if the particle is a gas particle then it has internal energy that must be added
#ifdef GASON
            if (Part[j+noffset].GetType()==GASTYPE)intE=Part[j+noffset].GetU();
#endif
#ifdef NOMASS
            Ti=(0.5*v2+intE)*opt.MassValue;
#else
            Ti=Part[j+noffset].GetMass()*(0.5*v2+intE);
#endif
            Part[j+noffset].SetPotential(Part[j+noffset].GetPotential()/Ti);
            pdata[i].T+=Ti;
            pdata[i].Ttyped[Part[j+noffset].GetType()]+=Ti;
            if(Part[j+noffset].GetPotential()<-1.0) {
                pdata[i].Efrac+=1.0;
                pdata[i].Efractyped[Part[j+noffset].GetType()]+=1.0;
            }
        }
        pdata[i].Efrac/=(Double_t)hp[i].AllNumberofParticles;
        for (j=0;j<NPARTTYPES;j++) if(hp[i].NumofType[j]>0) pdata[i].Efractyped[j]/=(Double_t)hp[i].AllNumofType[j];
        ngdone++;
    }
#ifdef USEOPENMP
}
#endif
    t1=MyGetTime()-t1;
    cout<<"Done "<<ngdone<<" small groups in "<<t1<<endl;
    t1=MyGetTime();
    //large groups use tree code, ids of the particles are overwritten to put particles back in original index order prior to gravity calculation
    //therefore must store for each large group the id values in pid values (which here)
    for (i=0;i<ngroup;i++) if (hp[i].AllNumberofParticles>=ompunbindnum){
        noffset=hp[i].noffset;
        if (opt.ipotcalc) Potential(opt,hp[i].NumberofParticles,&Part[noffset]);
        Tval=0;Potval=0;Efracval=0;
        for (j=0;j<omparraysize;j++) Tvaltyped[j]=Potvaltyped[j]=Efracvaltyped[j]=0.;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,v2,Ti,intE,poti,tid)
{
    #pragma omp for reduction(+:Tval,Efracval,Potval)
#endif
        for (j=0;j<hp[i].NumberofParticles;j++) {
#ifdef USEOPENMP
            tid=omp_get_thread_num();
#else
            tid=0;
#endif
            //reset id values
            //Part[noffset+j].SetID(Part[noffset+j].GetPID());
            v2=0.;for (int n=0;n<3;n++) v2+=pow(Part[j+noffset].GetVelocity(n),2.0);
            intE=0;
            //if the particle is a gas particle then it has internal energy that must be added
#ifdef GASON
            if (Part[j+noffset].GetType()==GASTYPE)intE=Part[j+noffset].GetU();
#endif
#ifdef NOMASS
            Tval+=Ti=(0.5*v2+intE)*opt.MassValue;
            Potval+=poti=Part[j+noffset].GetPotential()*mv2;
#else
            Tval+=Ti=(0.5*v2+intE)*Part[j+noffset].GetMass();
            Potval+=poti=Part[j+noffset].GetPotential();
#endif
            Part[j+noffset].SetPotential(poti/Ti);
            Tvaltyped[Part[j+noffset].GetType()+tid*NPARTTYPES]+=Ti;
            Potvaltyped[Part[j+noffset].GetType()+tid*NPARTTYPES]+=poti;
            if(Part[j+noffset].GetPotential()<-1.0) {
                Efracval+=1.0;
                Efracvaltyped[Part[j+noffset].GetType()+tid*NPARTTYPES]+=1.0;
            }
        }
#ifdef USEOPENMP
}
#endif
        pdata[i].T=Tval;pdata[i].Efrac=Efracval;pdata[i].Pot=Potval;
        pdata[i].Efrac/=(Double_t)hp[i].NumberofParticles;
        for (j=0;j<NPARTTYPES;j++) if(hp[i].NumofType[j]>0) {
            pdata[i].Ttyped[j]=pdata[i].Pottyped[j]=pdata[i].Efractyped[j]=0.0;
            for (k=0;k<nthreads;k++) {
                pdata[i].Ttyped[j]+=Tvaltyped[j+k*NPARTTYPES];
                pdata[i].Pottyped[j]+=Potvaltyped[j+k*NPARTTYPES];
                pdata[i].Efractyped[j]+=Efracvaltyped[j+k*NPARTTYPES];
            }
            pdata[i].Efractyped[j]/=(Double_t)hp[i].NumofType[j];
        }
    }
    t1=MyGetTime()-t1;
    cout<<"Done large groups in "<<t1<<endl;
}

///just calculates the potential energy of a particle
void GetPotentialEnergy(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp) {
    cout<<"Get Potential Energy"<<endl;
    Particle *Pval;
    Int_t i,j,k;
    Double_t eps2=opt.uinfo.eps*opt.uinfo.eps;
    Double_t r2,poti,pot;
    Double_t Potval;
    Double_t *Potvaltyped;
    Double_t mv2=opt.MassValue*opt.MassValue;
    Int_t noffset;
    Double_t t1;
    Int_t ngdone=0;
    int maxnthreads=1,nthreads=1,tid;
    Int_t omparraysize;
#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=nthreads=omp_get_num_threads();
    }
#endif
    omparraysize=nthreads*NPARTTYPES;
    Potvaltyped=new Double_t[omparraysize];

    t1=MyGetTime();
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,r2,poti,pot,noffset)
{
    #pragma omp for schedule(dynamic,1) reduction(+:ngdone) nowait
#endif

    for (i=0;i<ngroup;i++) if (hp[i].AllNumberofParticles<ompunbindnum && hp[i].AllNumberofParticles>0) {
        noffset=hp[i].noffset;
	pdata[i].Pot=0;
        for (j=0;j<NPARTTYPES;j++) pdata[i].Pottyped[j]=0.;
        for (j=0;j<hp[i].AllNumberofParticles;j++) {
            for (k=j+1;k<hp[i].AllNumberofParticles;k++) {
                r2=0.;for (int n=0;n<3;n++) r2+=pow(Part[j+noffset].GetPosition(n)-Part[k+noffset].GetPosition(n),2.0);
                r2+=eps2;
                r2=1.0/sqrt(r2);
#ifdef NOMASS
                pot=-opt.G*(Part[j+noffset].GetMass()*Part[k+noffset].GetMass())*r2*mv2;
#else
                pot=-opt.G*(Part[j+noffset].GetMass()*Part[k+noffset].GetMass())*r2;
#endif
                poti=Part[j+noffset].GetPotential()+pot;Part[j+noffset].SetPotential(poti);
                poti=Part[k+noffset].GetPotential()+pot;Part[k+noffset].SetPotential(poti);
            }
            pdata[i].Pot+=Part[j+noffset].GetPotential();
            pdata[i].Pottyped[Part[j+noffset].GetType()]+=Part[j+noffset].GetPotential();
        }
        ngdone++;
    }
#ifdef USEOPENMP
}
#endif
    t1=MyGetTime()-t1;
    cout<<"Done "<<ngdone<<" small groups in "<<t1<<endl;
    t1=MyGetTime();
    //large groups use tree code, ids of the particles are overwritten to put particles back in original index order prior to gravity calculation
    //therefore must store for each large group the id values in pid values (which here)
    for (i=0;i<ngroup;i++) if (hp[i].AllNumberofParticles>=ompunbindnum){
        noffset=hp[i].noffset;
        Potential(opt,hp[i].NumberofParticles,&Part[noffset]);
    }
    t1=MyGetTime()-t1;
    cout<<"Done large groups in "<<t1<<endl;
}

///Sort particles in according to their binding energy. Code very similar to that in STF/VELOCIraptor
//void SortAccordingtoBindingEnergy(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngroup, Int_t *pfof, Int_t *numingroup, Int_t *numingrouptype, PropData *pdata, Int_t ioffset)
void SortAccordingtoBindingEnergy(Options &opt, Particle *Part, Int_t ngroup, PropData **allpdata, HaloParticleData *hp, int icalcbindeflag, int itypesubsort,int ipotsort)
{
#ifndef USEMPI
    int ThisTask=0;
#endif
    cout<<ThisTask<<" Sort according to binding energy"<<endl;
    Int_t i,j,k;
    Int_t noffsettype;
    Double_t t1=MyGetTime();
    PropData *pdata=allpdata[DMTYPE];
    //then since dm centre of mass should have been loaded and corrected for periodicity
    //get binding energy of all particles in the group to the dm reference frame
    if (icalcbindeflag) {
        GetBindingEnergy(opt, Part, ngroup, pdata,hp);
        //store potential and kinetic info calculated by 
        for (i=0;i<ngroup;i++) {
            for (j=0;j<NPARTTYPES;j++) if (j!=DMTYPE&&hp[i].NumofType[j]>0) {
                allpdata[j][i].Efrac=pdata[i].Efrac;
                allpdata[j][i].T=pdata[i].T;
                allpdata[j][i].Pot=pdata[i].Pot;
                allpdata[j][i].Efractyped[j]=pdata[i].Efractyped[j];
                allpdata[j][i].Ttyped[j]=pdata[i].Ttyped[j];
                allpdata[j][i].Pottyped[j]=pdata[i].Pottyped[j];
            }
        }
    }
    if (itypesubsort)
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i)
{
    #pragma omp for nowait
#endif
        for (i=0;i<ngroup;i++)
            qsort(&Part[hp[i].noffset], hp[i].NumberofParticles, sizeof(Particle), TypeCompare);
#ifdef USEOPENMP
}
#endif
    //then since dm centre of mass should have been loaded and corrected for periodicity
    //get binding energy of all particles in the group to the dm reference frame
    //if (icalcbindeflag) GetBindingEnergy(opt, Part, ngroup, pdata,hp);
    cout<<"Sort particles by binding energy"<<endl;
    if (ipotsort) {
    if (itypesubsort) {
    cout<<"Sort by type and binding "<<endl; 
    //sort by type and energy
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,noffsettype)
{
    #pragma omp for nowait
#endif
    for (i=0;i<ngroup;i++) {
        noffsettype=0;
        for (j=0;j<NPARTTYPES;j++) {
            if (hp[i].NumofType[j]>0)
                qsort(&Part[hp[i].noffset+noffsettype], hp[i].NumofType[j], sizeof(Particle), PotCompare);
            noffsettype+=hp[i].NumofType[j];
        }
    }
#ifdef USEOPENMP
}
#endif
    }
    else {
    cout<<"sort all only by binding energy "<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i)
{
    #pragma omp for nowait
#endif
    for (i=0;i<ngroup;i++){
        qsort(&Part[hp[i].noffset], hp[i].NumberofParticles, sizeof(Particle), PotCompare);
    }
#ifdef USEOPENMP
}
#endif
    }
    }
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

///gets most bound particle of a specific type assuming particles have been ordered according to type then binding energy
void GetMostBoundParticle(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp)
{
    cout<<"Get Most bound particles "<<endl;
#ifndef USEMPI
    int ThisTask=0;
#endif
    Int_t i,j;
    Int_t noffsettype;
    Double_t t1=MyGetTime();

#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,noffsettype)
{
    #pragma omp for nowait
#endif
    for (i=0;i<ngroup;i++) if (hp[i].AllNumofType[pdata[i].ptype]>0){
        noffsettype=0;
        for (j=0;j<pdata[i].ptype;j++) noffsettype+=hp[i].AllNumofType[j];
        for (j=0;j<3;j++) {pdata[i].gpos[j]=Part[hp[i].noffset+noffsettype].GetPosition(j);pdata[i].gvel[j]=Part[hp[i].noffset+noffsettype].GetVelocity(j);}
    }
#ifdef USEOPENMP
}
#endif
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

///adjusts the Particle array so that all unbound particles are at the end, and the bound particle list is sorted
///by type. Adjusts the \ref HaloParticleData structure so that \ref NumberofParticles and \ref NumofType reflect the new bound particle lists
void GetBoundParticles(Options &opt, Particle *Part, Int_t ngroup, PropData **allpdata, HaloParticleData *hp, Double_t TVratio)
{
    cout<<"Get Bound Particles "<<endl;
#ifndef USEMPI
    int ThisTask=0;
#endif
    int maxnthreads=1,nthreads=1,tid;
    Int_t i,j,k,l;
    Int_t npt,*np,noffsettype,noffsettypeold;
    Double_t t1=MyGetTime();
    SortAccordingtoBindingEnergy(opt, Part, ngroup, allpdata, hp,0,1,1);

#ifdef USEOPENMP
#pragma omp parallel
    {
    if (omp_get_thread_num()==0) maxnthreads=nthreads=omp_get_num_threads();
    }
#endif
    np=new Int_t[NPARTTYPES*nthreads];
    cout<<"checking bound particles "<<endl;

#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,l,npt,noffsettype,noffsettypeold,tid)
{
#endif
#ifdef USEOPENMP
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (hp[i].AllNumberofParticles>0){
        npt=0;
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        for (k=0;k<NPARTTYPES;k++) np[k+tid*NPARTTYPES]=0;
        noffsettype=0;
        for (k=0;k<NPARTTYPES;k++) {
            if (k>0)noffsettype+=hp[i].AllNumofType[k-1];
            j=noffsettype+hp[i].noffset;
            l=0;
            while (Part[j].GetPotential() < TVratio && l<hp[i].NumofType[k]) {
                np[Part[j].GetType()+tid*NPARTTYPES]++;
                npt++;
                j++;l++;
            }
        }
        noffsettype=0;
        noffsettypeold=0;
        hp[i].NumberofParticles=npt;
        for (k=0;k<NPARTTYPES;k++) hp[i].NumofType[k]=np[k+tid*NPARTTYPES];
    }
#ifdef USEOPENMP
}
#endif
    delete[] np;

    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}
