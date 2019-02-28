/*! \file unbind.cxx
 *  \brief this file contains routines to check if groups are self-bound and if not unbind them as requried

    \todo Need to improve the gravity calculation (ie: for tree potential use something other than just monopole and also apply corrections if necessary).
    \todo Need to clean up unbind proceedure, ensure its mpi compatible and can be combined with a pglist output easily
 */

#include "stf.h"

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

///\name Remove unbound particles from a candidate group
//@{
/*!
    Interface for unbinding proceedure. Unbinding routine requires several arrays, such as numingroup, pglist,gPart,ids, etc
    This arrays may have been constructed prior to the unbinding call and so can be passed to the routine
    if this is called it uses Particle array then deletes it.
*/
int CheckUnboundGroups(Options opt, const Int_t nbodies, Particle *Part, Int_t &ngroup, Int_t *&pfof, Int_t *numingroup, Int_t **pglist, int ireorder, Int_t *groupflag){
    bool ningflag=false, pglistflag=false;
    int iflag;
    Int_t ng=ngroup;
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
    qsort(Part,nbodies,sizeof(Particle),IDCompare);
    noffset[0]=noffset[1]=0;
    for (Int_t i=2;i<=ngroup;i++) noffset[i]=noffset[i-1]+numingroup[i-1];
#else
    Particle **gPart=BuildPartList(ngroup, numingroup, pglist, Part);
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
    qsort(Part,nbodies,sizeof(Particle),PIDCompare);
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
    if (groupflag!=NULL) iflag=Unbind(opt, gPart, ngroup, numingroup,pfof,pglist,0);
    else iflag=Unbind(opt, gPart, ngroup, numingroup,pfof,pglist,ireorder);
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
            if (ireorder) ReorderGroupIDsAndArraybyValue(ng,ngroup,numingroup,pfof,pglist,groupvalue,groupflag);
            delete[] groupvalue;
        }
        else ReorderGroupIDs(ng,ngroup,numingroup,pfof,pglist);
    }
    for (Int_t i=1;i<=ng;i++) delete[] gPart[i];delete[] gPart;
#endif
    if (pglistflag) {for (Int_t i=1;i<=ng;i++) delete[] pglist[i];delete[] pglist;}
    if (ningflag) delete[] numingroup;

    if (opt.iverbose) cout<<ThisTask<<" Done. Number of groups remaining "<<ngroup<<endl;

    return iflag;
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
    int iunbindsizeflag;
    int maxnthreads,nthreads=1,l,n;
    Int_t i,j,k,ng=numgroups;
    Double_t maxE,totT,v2,r2,poti,Ti,eps2=opt.uinfo.eps*opt.uinfo.eps,mv2=opt.MassValue*opt.MassValue,Efrac;
    Double_t *gmass,*totV;
    PriorityQueue *pq;
    Int_t nEplus,pqsize, nEfrac;
    Int_t *nEplusid;
    int *Eplusflag;
    bool unbindcheck;
    Coordinate *cmvel;
    //for tree code potential calculation
    KDTree *tree;
    Int_t ncell,ntreecell,nleafcell;
    Int_t *start,*end;
    Double_t *cmtot,*cBmax,*cR2max, **r2val;
    Coordinate *cellcm;
    Node *root,**nodelist, **npomp;
    Int_t **marktreecell,**markleafcell;

    //used to determine potential based reference velocity frame
    Double_t potmin,menc;
    Int_t npot,ipotmin;
    Coordinate potpos;
    Int_t *storeval;

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    //note that it is possible that called as a library, velociraptor
    //does not need to calculate potentials itself
    //in that case do not calculate potentials but instead
    //copy relevant information
    cmvel   =new Coordinate[numgroups+1];
    gmass   =new Double_t[numgroups+1];
    totV    =new Double_t[numgroups+1];
    for (i=1;i<=numgroups;i++) {
        cmvel[i]=Coordinate(0.);
        gmass[i]=totV[i]=0.;
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
    if (opt.uinfo.icalculatepotential) {
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
#pragma omp parallel default(shared)  \
private(i,j,k,n,r2,poti)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=1;i<=numgroups;i++)
    {
        if (numingroup[i]<=UNBINDNUM) {
            for (j=0;j<numingroup[i];j++) {
                for (k=j+1;k<numingroup[i];k++) {
                    r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[i][j].GetPosition(n)-gPart[i][k].GetPosition(n),2.0);
                    r2+=eps2;
                    r2=1.0/sqrt(r2);
                    Double_t pot=-opt.G*(gPart[i][j].GetMass()*gPart[i][k].GetMass())*r2;
                    poti=gPart[i][j].GetPotential()+pot;gPart[i][j].SetPotential(poti);
                    poti=gPart[i][k].GetPotential()+pot;gPart[i][k].SetPotential(poti);
                }
            }
#ifdef NOMASS
        for (j=0;j<numingroup[i];j++) gPart[i][j].SetPotential(gPart[i][j].GetPotential()*mv2);
#endif
        for (j=0;j<numingroup[i];j++) totV[i]+=0.5*gPart[i][j].GetPotential();
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

    //now begin large group calculation
    marktreecell=new Int_t*[nthreads];
    markleafcell=new Int_t*[nthreads];
    r2val=new Double_t*[nthreads];
    npomp=new Node*[nthreads];
    //otherwise use tree tree gravity calculation
    //here openmp is per group since each group is large
    for (i=1;i<=numgroups;i++)
    {
        if (numingroup[i]>UNBINDNUM) {
            //to make this memory efficient really need just KDTree that uses Coordinates
            tree=new KDTree(gPart[i],numingroup[i],opt.uinfo.BucketSize,tree->TPHYS);

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
                    for (n=0;n<3;n++) cellcm[j][n]+=gPart[i][k].GetPosition(n)*gPart[i][k].GetMass();
                    cmtot[j]+=gPart[i][k].GetMass();
                }
                for (n=0;n<3;n++) cellcm[j][n]/=cmtot[j];
                Double_t xdiff,xdiff1;
                xdiff=(cellcm[j]-Coordinate(gPart[i][start[j]].GetPosition())).Length();
                for (k=start[j]+1;k<end[j];k++) {
                    xdiff1=(cellcm[j]-Coordinate(gPart[i][k].GetPosition())).Length();
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
private(j,k,l,n,ntreecell,nleafcell,r2,poti)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
            for (j=0;j<numingroup[i];j++) {
                int tid;
#ifdef USEOPENMP
                tid=omp_get_thread_num();
#else
                tid=0;
#endif
                npomp[tid]=tree->GetRoot();
                Coordinate xpos(gPart[i][j].GetPosition());
                nleafcell=ntreecell=0;
                MarkCell(npomp[tid],marktreecell[tid], markleafcell[tid],ntreecell,nleafcell,r2val[tid],opt.uinfo.BucketSize, cR2max, cellcm, cmtot, xpos, eps2);
                poti=0;
                for (k=0;k<ntreecell;k++) {
                    poti+=-gPart[i][j].GetMass()*r2val[tid][k];
                }
                for (k=0;k<nleafcell;k++) {
                    for (l=start[markleafcell[tid][k]];l<end[markleafcell[tid][k]];l++) {
                        if (j!=l) {
                            r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[i][j].GetPosition(n)-gPart[i][l].GetPosition(n),(Double_t)2.0);
                            r2+=eps2;
                            r2=1.0/sqrt(r2);
                            poti+=-(gPart[i][j].GetMass()*gPart[i][l].GetMass())*r2;
                        }
                    }
                }
                poti*=opt.G;
#ifdef NOMASS
                poti*=mv2;
#endif
                gPart[i][j].SetPotential(poti);
            }
#ifdef USEOPENMP
}
#endif
            for (j=0;j<numingroup[i];j++) totV[i]+=0.5*gPart[i][j].GetPotential();
            delete tree;
            delete[] start;
            delete[] end;
            delete[] cmtot;
            delete[] cBmax;
            delete[] cR2max;
            delete[] cellcm;
            delete[] nodelist;
            for (j=0;j<nthreads;j++) {delete[] marktreecell[j];delete[] markleafcell[j];delete[] r2val[j];}
        }
    }

    }//end of check whether we calculate potential

    //Now set the kinetic reference frame
    //if using standard frame, then using CMVEL of the entire structure
    if (opt.uinfo.cmvelreftype==CMVELREF) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,npot,menc,potpos,storeval)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=1;i<=numgroups;i++)
    {
        for (k=0;k<3;k++) potpos[k]=cmvel[i][k]=0;
        for (j=0;j<numingroup[i];j++) {
            gmass[i]+=gPart[i][j].GetMass();
            for (k=0;k<3;k++)
                potpos[k]+=gPart[i][j].GetPosition(k)*gPart[i][j].GetMass();
                //cmvel[i][k]+=gPart[i][j].GetVelocity(k)*gPart[i][j].GetMass();
        }
        for (k=0;k<3;k++)potpos[k]*=(1.0/gmass[i]);
        //for (k=0;k<3;k++)cmvel[i][k]*=(1.0/gmass[i]);
        //now sort by radius, first store original position in array
        storeval=new Int_t[numingroup[i]];
        for (j=0;j<numingroup[i];j++) {
            for (k=0;k<3;k++) gPart[i][j].SetPosition(k,gPart[i][j].GetPosition(k)-potpos[k]);
        }
        //sort by radius
        gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),RadCompare);
        //use central regions to define centre of mass velocity
        //determine how many particles to use
        npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
        npot=min(npot,numingroup[i]);
        cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
        for (j=0;j<npot;j++) {
            for (k=0;k<3;k++) cmvel[i][k]+=gPart[i][j].GetVelocity(k)*gPart[i][j].GetMass();
            menc+=gPart[i][j].GetMass();
        }
        for (j=0;j<3;j++) {cmvel[i][j]/=menc;}
        gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),IDCompare);
        for (j=0;j<numingroup[i];j++) {
            gPart[i][j].SetID(storeval[j]);
            for (k=0;k<3;k++) gPart[i][j].SetPosition(k,gPart[i][j].GetPosition(k)+potpos[k]);
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
private(i,j,k,npot,menc,potmin,ipotmin,potpos,storeval)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i=1;i<=numgroups;i++) {
            //determine how many particles to use
            npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
            npot=min(npot,numingroup[i]);

            storeval=new Int_t[numingroup[i]];
            for (j=0;j<numingroup[i];j++) {storeval[j]=gPart[i][j].GetID();gPart[i][j].SetID(j);}
            //determine position of minimum potential and by radius around this position
            potmin=gPart[i][0].GetPotential();ipotmin=0;
            for (j=1;j<numingroup[i];j++) if (gPart[i][j].GetPotential()<potmin) {potmin=gPart[i][j].GetPotential();ipotmin=j;}
            for (k=0;k<3;k++) potpos[k]=gPart[i][ipotmin].GetPosition(k);

            for (j=0;j<numingroup[i];j++) {
                for (k=0;k<3;k++) gPart[i][j].SetPosition(k,gPart[i][j].GetPosition(k)-potpos[k]);
            }
            gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),RadCompare);
            //now determine kinetic frame
            cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
            for (j=0;j<npot;j++) {
                for (k=0;k<3;k++) cmvel[i][k]+=gPart[i][j].GetVelocity(k)*gPart[i][j].GetMass();
                menc+=gPart[i][j].GetMass();
            }
            for (j=0;j<3;j++) {cmvel[i][j]/=menc;}
            gsl_heapsort(gPart[i],numingroup[i],sizeof(Particle),IDCompare);
            for (j=0;j<numingroup[i];j++) gPart[i][j].SetID(storeval[j]);
            delete[] storeval;
        }
#ifdef USEOPENMP
}
#endif
    }


    //now go through groups and begin unbinding by finding least bound particle
    //again for small groups multithread over groups
    //larger groups thread over particles in a group
    //for large groups, paralleize over particle, for small groups parallelize over groups
    //here energy data is stored in density
    for (i=1;i<=numgroups;i++) if (numingroup[i]>=ompunbindnum)
    {
        totT=0;
        Efrac=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,k,v2,Ti,unbindcheck)
{
    #pragma omp for reduction(+:totT,Efrac)
#endif
        for (j=0;j<numingroup[i];j++) {
            v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[i][j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
            Ti=0.5*gPart[i][j].GetMass()*v2*opt.MassValue;
#ifdef GASON
            Ti+=opt.MassValue*gPart[i][j].GetU();
#endif
#else
            Ti=0.5*gPart[i][j].GetMass()*v2;
#ifdef GASON
            Ti+=gPart[i][j].GetMass()*gPart[i][j].GetU();
#endif
#endif
            totT+=Ti;
            gPart[i][j].SetDensity(opt.uinfo.Eratio*Ti+gPart[i][j].GetPotential());
            Efrac+=(Ti+gPart[i][j].GetPotential()<0);
        }
#ifdef USEOPENMP
}
#endif
        Efrac/=(Double_t)numingroup[i];
        //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
        nEplus=0;
        pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+2);
        nEplusid=new Int_t[pqsize];
        Eplusflag=new int[numingroup[i]];
        maxE=gPart[i][0].GetDensity();
        for (j=1;j<numingroup[i];j++) if(maxE<gPart[i][j].GetDensity()) maxE=gPart[i][j].GetDensity();
        //check if bound;
        if (opt.uinfo.unbindtype==USYSANDPART) {
            if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
            else unbindcheck=false;
        }
        else if (opt.uinfo.unbindtype==UPART) {
            if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
            else unbindcheck=false;
        }
        //if need to unbind load largest energies so as to remove at most pqsize particles (roughly 1%) per removal loop
        if (unbindcheck) {
            for (j=0;j<numingroup[i];j++)Eplusflag[j]=0;
            pq=new PriorityQueue(pqsize);
            for (j=0;j<pqsize;j++) pq->Push(j,gPart[i][j].GetDensity());
            for (j=pqsize;j<numingroup[i];j++) if (gPart[i][j].GetDensity()>gPart[i][pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[i][j].GetDensity());}
            nEplus=0;
            //if just looking at particle then add to removal list till energy >0
            if (opt.uinfo.unbindtype==UPART) {
                for (j=0;j<pqsize;j++) {
                    if (gPart[i][pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
            else if (opt.uinfo.unbindtype==USYSANDPART) {
                nEfrac=0;
                if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                for (j=0;j<pqsize;j++) {
                    if (gPart[i][pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            delete pq;
        }

        while(unbindcheck)
        {
            iunbindflag+=1;
            //first correct for removal of all least bound particle
            double temp=1.0/gmass[i], temp2=0.;
            if (opt.uinfo.cmvelreftype==CMVELREF) {
                for (j=0;j<nEplus;j++) {
                    for (k=0;k<3;k++) cmvel[i][k]-=gPart[i][nEplusid[j]].GetVelocity(k)*gPart[i][nEplusid[j]].GetMass()*temp;
                    temp2+=gPart[i][nEplusid[j]].GetMass();
                }
                temp=gmass[i]/(gmass[i]-temp2);
                for (k=0;k<3;k++) cmvel[i][k]*=temp;
                gmass[i]-=temp2;
            }
            else {
                for (j=0;j<nEplus;j++) gmass[i]-=gPart[i][nEplusid[j]].GetMass();
            }
            //if ignore the background then adjust the potential energy of the particles
            //for large groups with many particles removed more computationally effective to simply
            //recalculate the potential energy after removing particles
            //for smaller number of particles removed, simply remove the contribution of this particle
            //from all others. The change in efficiency occurs at roughly nEplus>~log(numingroup[i]) particles. Here
            //we set the limit at 2*log(numingroup[i]) to account for overhead in producing tree and calculating new potential
            iunbindsizeflag=(nEplus<2.0*log((double)numingroup[i]));
            if (iunbindsizeflag) {
                if (opt.uinfo.bgpot==0) {
                    for (k=0;k<nEplus;k++) {
                        totV[i]-=0.5*gPart[i][nEplusid[k]].GetPotential();
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,r2,poti)
{
    #pragma omp for schedule(dynamic) nowait
#endif
                        for (j=0;j<numingroup[i];j++) {
                            if (j!=nEplusid[k]) {
                                r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[i][nEplusid[k]].GetPosition(n)-gPart[i][j].GetPosition(n),2.0);
                                r2+=eps2;
                                r2=1.0/sqrt(r2);
#ifdef NOMASS
                                poti=gPart[i][j].GetPotential()+opt.G*(gPart[i][nEplusid[k]].GetMass()*gPart[i][j].GetMass())*r2*mv2;
#else
                                poti=gPart[i][j].GetPotential()+opt.G*(gPart[i][nEplusid[k]].GetMass()*gPart[i][j].GetMass())*r2;
#endif
                                gPart[i][j].SetPotential(poti);
                            }
                        }
#ifdef USEOPENMP
}
#endif
                    }
                }
            }
            else {
                if (opt.uinfo.bgpot==0) for (k=0;k<nEplus;k++) totV[i]-=0.5*gPart[i][nEplusid[k]].GetPotential();
            }
            //remove particles with positive energy
            for (j=0;j<nEplus;j++) pfof[pglist[i][nEplusid[j]]]=0;
            k=numingroup[i]-1;
            for (j=0;j<nEplus;j++) if (nEplusid[j]<numingroup[i]-nEplus) {
                while(Eplusflag[k]==1)k--;
                pglist[i][nEplusid[j]]=pglist[i][k];
                gPart[i][nEplusid[j]]=gPart[i][k];
                Eplusflag[nEplusid[j]]=0;
                k--;
            }
            numingroup[i]-=nEplus;
            if (iunbindsizeflag && opt.uinfo.bgpot==0) Potential(opt, numingroup[i], gPart[i]);
            //if number of particles remove with positive energy is near to the number allowed to be removed
            //must recalculate kinetic energies and check if maxE>0
            //otherwise, end unbinding.
            if (nEplus>=0.1*pqsize+0.5) {

            //recalculate kinetic energies since cmvel has changed
            totT=0.;
            Efrac=0.;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,k,v2,Ti,unbindcheck)
{
    #pragma omp for reduction(+:totT,Efrac)
#endif
        for (j=0;j<numingroup[i];j++) {
            v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[i][j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
            Ti=0.5*gPart[i][j].GetMass()*v2*opt.MassValue;
#ifdef GASON
            Ti+=opt.MassValue*gPart[i][j].GetU();
#endif
#else
            Ti=0.5*gPart[i][j].GetMass()*v2;
#ifdef GASON
            Ti+=gPart[i][j].GetMass()*gPart[i][j].GetU();
#endif
#endif
            totT+=Ti;
            gPart[i][j].SetDensity(opt.uinfo.Eratio*Ti+gPart[i][j].GetPotential());
            Efrac+=(Ti+gPart[i][j].GetPotential()<0);
        }
#ifdef USEOPENMP
}
#endif
            Efrac/=numingroup[i];
            //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
            maxE=gPart[i][0].GetDensity();
            for (j=1;j<numingroup[i];j++) if(maxE<gPart[i][j].GetDensity()) maxE=gPart[i][j].GetDensity();
            pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+1);
            if (opt.uinfo.unbindtype==USYSANDPART) {
                if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
                else unbindcheck=false;
            }
            else if (opt.uinfo.unbindtype==UPART) {
                if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
                else unbindcheck=false;
            }
            if (unbindcheck) {
                pq=new PriorityQueue(pqsize);
                nEplus=0;
                for (j=0;j<pqsize;j++) pq->Push(j,gPart[i][j].GetDensity());
                for (j=pqsize;j<numingroup[i];j++) if (gPart[i][j].GetDensity()>gPart[i][pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[i][j].GetDensity());}
                //if just looking at particle then add to removal list till energy >0
                if (opt.uinfo.unbindtype==UPART) {
                    for (j=0;j<pqsize;j++) {
                        if (gPart[i][pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
                else if (opt.uinfo.unbindtype==USYSANDPART) {
                    nEfrac=0;
                    if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                    for (j=0;j<pqsize;j++) {
                        if (gPart[i][pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                delete pq;
            }
        }
        else unbindcheck=false;
        }
        //if group too small remove entirely
        if (numingroup[i]<opt.MinSize) {
            for (j=0;j<numingroup[i];j++) pfof[pglist[i][j]]=0;
            numingroup[i]=0;
            Efrac=0;
            iunbindflag++;
        }
        delete[] nEplusid;
        delete[] Eplusflag;
    }

    //now for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,n,maxE,pq,pqsize,nEplus,nEplusid,Eplusflag,totT,v2,Ti,unbindcheck,Efrac,nEfrac,r2,poti)
{
    #pragma omp for schedule(dynamic) nowait reduction(+:iunbindflag)
#endif
    for (i=1;i<=numgroups;i++) if (numingroup[i]<=ompunbindnum && numingroup[i]>0)
    {
        totT=0;
        maxE=-MAXVALUE;
        nEplus=0;
        Efrac=0.;
        for (j=0;j<numingroup[i];j++) {
            v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[i][j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
            Ti=0.5*gPart[i][j].GetMass()*v2*opt.MassValue;
#ifdef GASON
            Ti+=opt.MassValue*gPart[i][j].GetU();
#endif
#else
            Ti=0.5*gPart[i][j].GetMass()*v2;
#ifdef GASON
            Ti+=gPart[i][j].GetMass()*gPart[i][j].GetU();
#endif
#endif
            totT+=Ti;
            gPart[i][j].SetDensity(opt.uinfo.Eratio*Ti+gPart[i][j].GetPotential());
            Efrac+=(Ti+gPart[i][j].GetPotential()<0);
            if(maxE<gPart[i][j].GetDensity()) maxE=gPart[i][j].GetDensity();
        }
        Efrac/=(Double_t)numingroup[i];
        //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
        pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+2);
        Eplusflag=new int[numingroup[i]];
        nEplusid=new Int_t[pqsize];
        if (opt.uinfo.unbindtype==USYSANDPART) {
            if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
            else unbindcheck=false;
        }
        else if (opt.uinfo.unbindtype==UPART) {
            if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
            else unbindcheck=false;
        }
        if (unbindcheck) {
            for (j=0;j<numingroup[i];j++)Eplusflag[j]=0;
            pq=new PriorityQueue(pqsize);
            for (j=0;j<pqsize;j++) pq->Push(j,gPart[i][j].GetDensity());
            for (j=pqsize;j<numingroup[i];j++) if (gPart[i][j].GetDensity()>gPart[i][pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[i][j].GetDensity());}
            nEplus=0;
            //if just looking at particle then add to removal list till energy >0
            if (opt.uinfo.unbindtype==UPART) {
                for (j=0;j<pqsize;j++) {
                    if (gPart[i][pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
            else if (opt.uinfo.unbindtype==USYSANDPART) {
                nEfrac=0;
                if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                for (j=0;j<pqsize;j++) {
                    if (gPart[i][pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            delete pq;
        }
        while(unbindcheck)
        {
            iunbindflag++;
            //first correct for removal of all least bound particle
            double temp=1.0/gmass[i], temp2=0.;
            if (opt.uinfo.cmvelreftype==CMVELREF) {
                for (j=0;j<nEplus;j++) {
                    for (k=0;k<3;k++) cmvel[i][k]-=gPart[i][nEplusid[j]].GetVelocity(k)*gPart[i][nEplusid[j]].GetMass()*temp;
                    temp2+=gPart[i][nEplusid[j]].GetMass();
                }
                temp=gmass[i]/(gmass[i]-temp2);
                for (k=0;k<3;k++) cmvel[i][k]*=temp;
                gmass[i]-=temp2;
            }
            else {
                for (j=0;j<nEplus;j++) gmass[i]-=gPart[i][nEplusid[j]].GetMass();
            }
            //if ignore the background then adjust the potential energy of the particles
            if (opt.uinfo.bgpot==0) {
                for (k=0;k<nEplus;k++) {
                    totV[i]-=0.5*gPart[i][nEplusid[k]].GetPotential();
                    for (j=0;j<numingroup[i];j++) {
                        if (j!=nEplusid[k]) {
                            r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[i][nEplusid[k]].GetPosition(n)-gPart[i][j].GetPosition(n),2.0);
                            r2+=eps2;
                            r2=1.0/sqrt(r2);
#ifdef NOMASS
                            poti=gPart[i][j].GetPotential()+opt.G*(gPart[i][nEplusid[k]].GetMass()*gPart[i][j].GetMass())*r2*mv2;
#else
                            poti=gPart[i][j].GetPotential()+opt.G*(gPart[i][nEplusid[k]].GetMass()*gPart[i][j].GetMass())*r2;
#endif
                            gPart[i][j].SetPotential(poti);
                        }
                    }
                }
            }
            //remove particles with positive energy
            for (j=0;j<nEplus;j++) pfof[pglist[i][nEplusid[j]]]=0;
            k=numingroup[i]-1;
            for (j=0;j<nEplus;j++) if (nEplusid[j]<numingroup[i]-nEplus) {
                while(Eplusflag[k]==1)k--;
                pglist[i][nEplusid[j]]=pglist[i][k];
                gPart[i][nEplusid[j]]=gPart[i][k];
                Eplusflag[nEplusid[j]]=0;
                k--;
            }
            numingroup[i]-=nEplus;

            if (nEplus>=0.1*pqsize+0.5) {

            //recalculate kinetic energies since cmvel has changed
            totT=0;
            maxE=-MAXVALUE;
            for (j=0;j<numingroup[i];j++) {
                v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[i][j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
                Ti=0.5*gPart[i][j].GetMass()*v2*opt.MassValue;
#ifdef GASON
                Ti+=opt.MassValue*gPart[i][j].GetU();
#endif
#else
                Ti=0.5*gPart[i][j].GetMass()*v2;
#ifdef GASON
                Ti+=gPart[i][j].GetMass()*gPart[i][j].GetU();
#endif
#endif
                totT+=Ti;
                gPart[i][j].SetDensity(opt.uinfo.Eratio*Ti+gPart[i][j].GetPotential());
                Efrac+=(Ti+gPart[i][j].GetPotential()<0);
                    if(maxE<gPart[i][j].GetDensity()) maxE=gPart[i][j].GetDensity();
            }
            Efrac/=(Double_t)numingroup[i];
            if (opt.uinfo.unbindtype==USYSANDPART) {
                if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
                else unbindcheck=false;
            }
            else if (opt.uinfo.unbindtype==UPART) {
                if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
                else unbindcheck=false;
            }
            if (unbindcheck) {
                //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
                pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+1);
                pq=new PriorityQueue(pqsize);
                for (j=0;j<pqsize;j++) pq->Push(j,gPart[i][j].GetDensity());
                for (j=pqsize;j<numingroup[i];j++) if (gPart[i][j].GetDensity()>gPart[i][pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[i][j].GetDensity());}
                nEplus=0;
                //if just looking at particle then add to removal list till energy >0
                if (opt.uinfo.unbindtype==UPART) {
                    for (j=0;j<pqsize;j++) {
                        if (gPart[i][pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
                else if (opt.uinfo.unbindtype==USYSANDPART) {
                    nEfrac=0;
                    if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                    for (j=0;j<pqsize;j++) {
                        if (gPart[i][pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                delete pq;
            }
        }
        else unbindcheck=false;
        }
        //if group too small remove entirely
        if (numingroup[i]<opt.MinSize) {
            iunbindflag++;
            for (j=0;j<numingroup[i];j++) pfof[pglist[i][j]]=0;
            numingroup[i]=0;
            Efrac=0;
        }
        delete[] nEplusid;
        delete[] Eplusflag;
    }
#ifdef USEOPENMP
}
#endif

    for (i=1;i<=numgroups;i++) if (numingroup[i]==0) ng--;
    if (ireorder==1 && iunbindflag&&ng>0) ReorderGroupIDs(numgroups,ng,numingroup,pfof,pglist);
    delete[] cmvel;
    numgroups=ng;
    //return if any unbinding done indicating groups have been reordered
    if (iunbindflag) return 1;
    else return 0;
}

///Similar to unbind algorithm but assumes particles are ordered. saves memory but more computations
int Unbind(Options &opt, Particle *&gPart, Int_t &numgroups, Int_t *&numingroup, Int_t *&noffset, Int_t *&pfof)
{
    //flag which is changed if any groups are altered as groups may need to be reordered.
    int iunbindflag=0;
    //flag used to determine what style of update to the potential is done for larger groups as
    //if the amount of particles removed is large enough for large groups, it is more efficient to
    //recalculate the entire potential using a Tree code than it is removing the contribution of each removed particle from
    //all other particles
    int iunbindsizeflag;
    int maxnthreads,nthreads=1,l,n;
    Int_t i,j,k,ng=numgroups;
    Double_t maxE,totT,v2,r2,poti,Ti,eps2=opt.uinfo.eps*opt.uinfo.eps,mv2=opt.MassValue*opt.MassValue,Efrac;
    Double_t *gmass,*totV;
    PriorityQueue *pq;
    Int_t nEplus,pqsize,nEfrac;
    Int_t *nEplusid;
    int *Eplusflag;
    bool unbindcheck;
    Coordinate *cmvel;
    Particle Ptemp;

    //for tree code potential calculation
    KDTree *tree;
    Int_t ncell,ntreecell,nleafcell;
    Int_t *start,*end;
    Double_t *cmtot,*cBmax,*cR2max, **r2val;
    Coordinate *cellcm;
    Node *root,**nodelist, **npomp;
    Int_t **marktreecell,**markleafcell;

    //used to determine potential based reference velocity frame
    Double_t potmin,menc;
    Int_t npot,ipotmin;
    Coordinate potpos;
    Int_t *storeval;

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    cmvel   =new Coordinate[numgroups+1];
    gmass   =new Double_t[numgroups+1];
    totV    =new Double_t[numgroups+1];
    for (i=1;i<=numgroups;i++) {
        cmvel[i]=Coordinate(0.);
        gmass[i]=totV[i]=0.;
        if (opt.uinfo.icalculatepotential) {
            for (j=0;j<numingroup[i];j++) gPart[noffset[i]+j].SetPotential(0);
        }
        #ifdef SWIFTINTERFACE
        else {
            for (j=0;j<numingroup[i];j++) gPart[noffset[i]+j].SetPotential(gPart[noffset[i]+j].GetGravityPotential());
        }
        #endif
    }

    //if calculate potential
    if (opt.uinfo.icalculatepotential) {
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
#pragma omp parallel default(shared)  \
private(i,j,k,n,r2,poti)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=1;i<=numgroups;i++)
    {
        if (numingroup[i]<=UNBINDNUM) {
            for (j=0;j<numingroup[i];j++) {
                for (k=j+1;k<numingroup[i];k++) {
                    r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[noffset[i]+j].GetPosition(n)-gPart[noffset[i]+k].GetPosition(n),2.0);
                    r2+=eps2;
                    r2=1.0/sqrt(r2);
                    Double_t pot=-opt.G*(gPart[noffset[i]+j].GetMass()*gPart[noffset[i]+k].GetMass())*r2;
                    poti=gPart[noffset[i]+j].GetPotential()+pot;gPart[noffset[i]+j].SetPotential(poti);
                    poti=gPart[noffset[i]+k].GetPotential()+pot;gPart[noffset[i]+k].SetPotential(poti);
                }
            }
#ifdef NOMASS
        for (j=0;j<numingroup[i];j++) gPart[noffset[i]+j].SetPotential(gPart[noffset[i]+j].GetPotential()*mv2);
#endif
        for (j=0;j<numingroup[i];j++) totV[i]+=0.5*gPart[noffset[i]+j].GetPotential();
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

    //now begin large group calculation
    marktreecell=new Int_t*[nthreads];
    markleafcell=new Int_t*[nthreads];
    r2val=new Double_t*[nthreads];
    npomp=new Node*[nthreads];
    //otherwise use tree tree gravity calculation
    //here openmp is per group since each group is large
    for (i=1;i<=numgroups;i++)
    {
        if (numingroup[i]>UNBINDNUM) {
            //to make this memory efficient really need just KDTree that uses Coordinates
            tree=new KDTree(&gPart[noffset[i]],numingroup[i],opt.uinfo.BucketSize,tree->TPHYS);

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
                    for (n=0;n<3;n++) cellcm[j][n]+=gPart[noffset[i]+k].GetPosition(n)*gPart[noffset[i]+k].GetMass();
                    cmtot[j]+=gPart[noffset[i]+k].GetMass();
                }
                for (n=0;n<3;n++) cellcm[j][n]/=cmtot[j];
                Double_t xdiff,xdiff1;
                xdiff=(cellcm[j]-Coordinate(gPart[noffset[i]+start[j]].GetPosition())).Length();
                for (k=start[j]+1;k<end[j];k++) {
                    xdiff1=(cellcm[j]-Coordinate(gPart[noffset[i]+k].GetPosition())).Length();
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
private(j,k,l,n,ntreecell,nleafcell,r2,poti)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
            for (j=0;j<numingroup[i];j++) {
                int tid;
#ifdef USEOPENMP
                tid=omp_get_thread_num();
#else
                tid=0;
#endif
                npomp[tid]=tree->GetRoot();
                Coordinate xpos(gPart[noffset[i]+j].GetPosition());
                nleafcell=ntreecell=0;
                MarkCell(npomp[tid],marktreecell[tid], markleafcell[tid],ntreecell,nleafcell,r2val[tid],opt.uinfo.BucketSize, cR2max, cellcm, cmtot, xpos, eps2);
                poti=0;
                for (k=0;k<ntreecell;k++) {
                    poti+=-gPart[noffset[i]+j].GetMass()*r2val[tid][k];
                }
                for (k=0;k<nleafcell;k++) {
                    for (l=start[markleafcell[tid][k]];l<end[markleafcell[tid][k]];l++) {
                        if (j!=l) {
                            r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[noffset[i]+j].GetPosition(n)-gPart[noffset[i]+l].GetPosition(n),(Double_t)2.0);
                            r2+=eps2;
                            r2=1.0/sqrt(r2);
                            poti+=-(gPart[noffset[i]+j].GetMass()*gPart[noffset[i]+l].GetMass())*r2;
                        }
                    }
                }
                poti*=opt.G;
#ifdef NOMASS
                poti*=mv2;
#endif
                gPart[noffset[i]+j].SetPotential(poti);
            }
#ifdef USEOPENMP
}
#endif
            for (j=0;j<numingroup[i];j++) totV[i]+=0.5*gPart[noffset[i]+j].GetPotential();
            delete tree;
            delete[] start;
            delete[] end;
            delete[] cmtot;
            delete[] cBmax;
            delete[] cR2max;
            delete[] cellcm;
            delete[] nodelist;
            for (j=0;j<nthreads;j++) {delete[] marktreecell[j];delete[] markleafcell[j];delete[] r2val[j];}
        }
    }
    }//end of if calculate potential

    //Now set the kinetic reference frame
    //if using standard frame, then using CMVEL of the entire structure
    if (opt.uinfo.cmvelreftype==CMVELREF) {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,npot,menc,potpos,storeval)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=1;i<=numgroups;i++)
    {
        for (j=0;j<numingroup[i];j++) {
            gmass[i]+=gPart[noffset[i]+j].GetMass();
            for (k=0;k<3;k++)
                cmvel[i][k]+=gPart[noffset[i]+j].GetVelocity(k)*gPart[noffset[i]+j].GetMass();
        }
        for (k=0;k<3;k++)cmvel[i][k]*=(1.0/gmass[i]);
        for (k=0;k<3;k++) potpos[k]=cmvel[i][k]=0;
        for (j=0;j<numingroup[i];j++) {
            gmass[i]+=gPart[noffset[i]+j].GetMass();
            for (k=0;k<3;k++)
                potpos[k]+=gPart[noffset[i]+j].GetPosition(k)*gPart[noffset[i]+j].GetMass();
        }
        for (k=0;k<3;k++)potpos[k]*=(1.0/gmass[i]);
        //for (k=0;k<3;k++)cmvel[i][k]*=(1.0/gmass[i]);
        //now sort by radius, first store original position in array
        storeval=new Int_t[numingroup[i]];
        for (j=0;j<numingroup[i];j++) {
            for (k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)-potpos[k]);
        }
        //sort by radius
        gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),RadCompare);
        //use central regions to define centre of mass velocity
        //determine how many particles to use
        npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
        npot=min(npot,numingroup[i]);
        cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
        for (j=0;j<npot;j++) {
            for (k=0;k<3;k++) cmvel[i][k]+=gPart[noffset[i]+j].GetVelocity(k)*gPart[noffset[i]+j].GetMass();
            menc+=gPart[noffset[i]+j].GetMass();
        }
        for (j=0;j<3;j++) {cmvel[i][j]/=menc;}
        gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),IDCompare);
        for (j=0;j<numingroup[i];j++) {
            gPart[noffset[i]+j].SetID(storeval[j]);
            for (k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)+potpos[k]);
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
private(i,j,k,npot,menc,potmin,ipotmin,potpos,storeval)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i=1;i<=numgroups;i++) {
            //must store original order
            storeval=new Int_t[numingroup[i]];
            for (j=0;j<numingroup[i];j++) {storeval[j]=gPart[noffset[i]+j].GetID();gPart[noffset[i]+j].SetID(j);}
            //determine how many particles to use
            npot=max(opt.uinfo.Npotref,Int_t(opt.uinfo.fracpotref*numingroup[i]));
            //determine position of minimum potential and by radius around this position
            potmin=gPart[noffset[i]+0].GetPotential();ipotmin=0;
            for (j=1;j<numingroup[i];j++) if (gPart[noffset[i]+j].GetPotential()<potmin) {potmin=gPart[noffset[i]+j].GetPotential();ipotmin=j;}
            for (k=0;k<3;k++) potpos[k]=gPart[noffset[i]+ipotmin].GetPosition(k);
            for (j=0;j<numingroup[i];j++) {
                for (k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)-potpos[k]);
            }
            gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),RadCompare);
            //now determine kinetic frame
            cmvel[i][0]=cmvel[i][1]=cmvel[i][2]=menc=0.;
            for (j=0;j<npot;j++) {
                for (k=0;k<3;k++) cmvel[i][k]+=gPart[noffset[i]+j].GetVelocity(k)*gPart[noffset[i]+j].GetMass();
                menc+=gPart[noffset[i]+j].GetMass();
            }
            for (j=0;j<3;j++) {cmvel[i][j]/=menc;}
            gsl_heapsort(&gPart[noffset[i]],numingroup[i],sizeof(Particle),IDCompare);
            for (j=0;j<numingroup[i];j++) {
                gPart[noffset[i]+j].SetID(storeval[j]);
                for (k=0;k<3;k++) gPart[noffset[i]+j].SetPosition(k,gPart[noffset[i]+j].GetPosition(k)+potpos[k]);
            }
            delete[] storeval;
        }
#ifdef USEOPENMP
}
#endif
    }

    //now go through groups and begin unbinding by finding least bound particle
    //again for small groups multithread over groups
    //larger groups thread over particles in a group
    //for large groups, paralleize over particle, for small groups parallelize over groups
    //here energy data is stored in density
    for (i=1;i<=numgroups;i++) if (numingroup[i]>=ompunbindnum)
    {
        totT=0;
        Efrac=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,k,v2,Ti,unbindcheck)
{
    #pragma omp for reduction(+:totT,Efrac)
#endif
        for (j=0;j<numingroup[i];j++) {
            v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[noffset[i]+j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
            Ti=0.5*gPart[noffset[i]+j].GetMass()*v2*opt.MassValue;
#ifdef GASON
            Ti+=opt.MassValue*gPart[noffset[i]+j].GetU();
#endif
#else
            Ti=0.5*gPart[noffset[i]+j].GetMass()*v2;
#ifdef GASON
            Ti+=gPart[noffset[i]+j].GetMass()*gPart[noffset[i]+j].GetU();
#endif
#endif
            totT+=Ti;
            gPart[noffset[i]+j].SetDensity(opt.uinfo.Eratio*Ti+gPart[noffset[i]+j].GetPotential());
            Efrac+=(Ti+gPart[noffset[i]+j].GetPotential()<0);
        }
#ifdef USEOPENMP
}
#endif
        Efrac/=(Double_t)numingroup[i];
        //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
        nEplus=0;
        pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+2);
        nEplusid=new Int_t[pqsize];
        Eplusflag=new int[numingroup[i]];
        maxE=gPart[noffset[i]+0].GetDensity();
        for (j=1;j<numingroup[i];j++) if(maxE<gPart[noffset[i]+j].GetDensity()) maxE=gPart[noffset[i]+j].GetDensity();
        //if largest energy is positive load largest energies so as to remove at most pqsize particles (roughly 1%) per removal loop
        if (opt.uinfo.unbindtype==USYSANDPART) {
            if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
            else unbindcheck=false;
        }
        else if (opt.uinfo.unbindtype==UPART) {
            if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
            else unbindcheck=false;
        }
        if (unbindcheck) {
            for (j=0;j<numingroup[i];j++)Eplusflag[j]=0;
            pq=new PriorityQueue(pqsize);
            for (j=0;j<pqsize;j++) pq->Push(j,gPart[noffset[i]+j].GetDensity());
            for (j=pqsize;j<numingroup[i];j++) if (gPart[noffset[i]+j].GetDensity()>gPart[noffset[i]+pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[noffset[i]+j].GetDensity());}
            nEplus=0;
            //if just looking at particle then add to removal list till energy >0
            if (opt.uinfo.unbindtype==UPART) {
                for (j=0;j<pqsize;j++) {
                    if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
            else if (opt.uinfo.unbindtype==USYSANDPART) {
                nEfrac=0;
                if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                for (j=0;j<pqsize;j++) {
                    if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            delete pq;
        }
        while(unbindcheck)
        {
            iunbindflag++;
            //first correct for removal of all least bound particle
            double temp=1.0/gmass[i], temp2=0.;
            if (opt.uinfo.cmvelreftype==CMVELREF) {
                for (j=0;j<nEplus;j++) {
                    for (k=0;k<3;k++) cmvel[i][k]-=gPart[noffset[i]+nEplusid[j]].GetVelocity(k)*gPart[noffset[i]+nEplusid[j]].GetMass()*temp;
                    temp2+=gPart[noffset[i]+nEplusid[j]].GetMass();
                }
                temp=gmass[i]/(gmass[i]-temp2);
                for (k=0;k<3;k++) cmvel[i][k]*=temp;
                gmass[i]-=temp2;
            }
            else {
                for (j=0;j<nEplus;j++) gmass[i]-=gPart[noffset[i]+nEplusid[j]].GetMass();
            }
            //if ignore the background then adjust the potential energy of the particles
            //for large groups with many particles removed more computationally effective to simply
            //recalculate the potential energy after removing particles
            //for smaller number of particles removed, simply remove the contribution of this particle
            //from all others. The change in efficiency occurs at roughly nEplus>~log(numingroup[i]) particles. Here
            //we set the limit at 2*log(numingroup[i]) to account for overhead in producing tree and calculating new potential
            iunbindsizeflag=(nEplus<2.0*log((double)numingroup[i]));
            if (iunbindsizeflag) {
                if (opt.uinfo.bgpot==0) {
                    for (k=0;k<nEplus;k++) {
                        totV[i]-=0.5*gPart[noffset[i]+nEplusid[k]].GetPotential();
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,r2,poti)
{
    #pragma omp for schedule(dynamic) nowait
#endif
                        for (j=0;j<numingroup[i];j++) {
                            if (j!=nEplusid[k]) {
                                r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[noffset[i]+nEplusid[k]].GetPosition(n)-gPart[noffset[i]+j].GetPosition(n),2.0);
                                r2+=eps2;
                                r2=1.0/sqrt(r2);
#ifdef NOMASS
                                poti=gPart[noffset[i]+j].GetPotential()+opt.G*(gPart[noffset[i]+nEplusid[k]].GetMass()*gPart[noffset[i]+j].GetMass())*r2*mv2;
#else
                                poti=gPart[noffset[i]+j].GetPotential()+opt.G*(gPart[noffset[i]+nEplusid[k]].GetMass()*gPart[noffset[i]+j].GetMass())*r2;
#endif
                                gPart[noffset[i]+j].SetPotential(poti);
                            }
                        }
#ifdef USEOPENMP
}
#endif
                    }
                }
            }
            else {
                if (opt.uinfo.bgpot==0) for (k=0;k<nEplus;k++) totV[i]-=0.5*gPart[noffset[i]+nEplusid[k]].GetPotential();
            }
            //remove particles with positive energy
            for (j=0;j<nEplus;j++) pfof[gPart[noffset[i]+nEplusid[j]].GetPID()]=0;
            k=numingroup[i]-1;
            for (j=0;j<nEplus;j++) if (nEplusid[j]<numingroup[i]-nEplus) {
                while(Eplusflag[k]==1)k--;
                Ptemp=gPart[noffset[i]+nEplusid[j]];
                gPart[noffset[i]+nEplusid[j]]=gPart[noffset[i]+k];
                gPart[noffset[i]+k]=Ptemp;
                Eplusflag[nEplusid[j]]=0;
                k--;
            }
            numingroup[i]-=nEplus;
            if (iunbindsizeflag && opt.uinfo.bgpot==0) Potential(opt, numingroup[i], &gPart[noffset[i]]);
            //if number of particles remove with positive energy is near to the number allowed to be removed
            //must recalculate kinetic energies and check if maxE>0
            //otherwise, end unbinding.
            if (nEplus>=0.1*pqsize+0.5) {

            //recalculate kinetic energies since cmvel has changed
            totT=0.;
            Efrac=0.;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,k,v2,Ti,unbindcheck)
{
    #pragma omp for reduction(+:totT,Efrac)
#endif
        for (j=0;j<numingroup[i];j++) {
            v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[noffset[i]+j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
            Ti=0.5*gPart[noffset[i]+j].GetMass()*v2*opt.MassValue;
#ifdef GASON
            Ti+=opt.MassValue*gPart[noffset[i]+j].GetU();
#endif
#else
            Ti=0.5*gPart[noffset[i]+j].GetMass()*v2;
#ifdef GASON
            Ti+=gPart[noffset[i]+j].GetMass()*gPart[noffset[i]+j].GetU();
#endif
#endif
            totT+=Ti;
            gPart[noffset[i]+j].SetDensity(opt.uinfo.Eratio*Ti+gPart[noffset[i]+j].GetPotential());
            Efrac+=(Ti+gPart[noffset[i]+j].GetPotential()<0);
        }
#ifdef USEOPENMP
}
#endif
            //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
            maxE=gPart[noffset[i]+0].GetDensity();
            for (j=1;j<numingroup[i];j++) if(maxE<gPart[noffset[i]+j].GetDensity()) maxE=gPart[noffset[i]+j].GetDensity();
            pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+1);
            if (opt.uinfo.unbindtype==USYSANDPART) {
                if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
                else unbindcheck=false;
            }
            else if (opt.uinfo.unbindtype==UPART) {
                if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
                else unbindcheck=false;
            }
            if (unbindcheck) {
                pq=new PriorityQueue(pqsize);
                nEplus=0;
                for (j=0;j<pqsize;j++) pq->Push(j,gPart[noffset[i]+j].GetDensity());
                for (j=pqsize;j<numingroup[i];j++) if (gPart[noffset[i]+j].GetDensity()>gPart[noffset[i]+pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[noffset[i]+j].GetDensity());}
                //if just looking at particle then add to removal list till energy >0
                if (opt.uinfo.unbindtype==UPART) {
                    for (j=0;j<pqsize;j++) {
                        if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
                else if (opt.uinfo.unbindtype==USYSANDPART) {
                    nEfrac=0;
                    if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                    for (j=0;j<pqsize;j++) {
                        if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                delete pq;
            }
            }
            else unbindcheck=false;
        }
        //if group too small remove entirely
        if (numingroup[i]<opt.MinSize) {
            iunbindflag++;
            for (j=0;j<numingroup[i];j++) pfof[gPart[noffset[i]+j].GetPID()]=0;
            numingroup[i]=0;
        }
        delete[] nEplusid;
        delete[] Eplusflag;
    }
    //now for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,maxE,pq,pqsize,nEplus,nEplusid,Eplusflag,totT,v2,Ti,unbindcheck,Efrac,Ptemp,nEfrac)
{
    #pragma omp for schedule(dynamic) nowait reduction(+:iunbindflag)
#endif
    for (i=1;i<=numgroups;i++) if (numingroup[i]<=ompunbindnum)
    {
        totT=0;
        maxE=-MAXVALUE;
        nEplus=0;
        Efrac=0.;
        for (j=0;j<numingroup[i];j++) {
            v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[noffset[i]+j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
            Ti=0.5*gPart[noffset[i]+j].GetMass()*v2*opt.MassValue;
#ifdef GASON
            Ti+=opt.MassValue*gPart[noffset[i]+j].GetU();
#endif
#else
            Ti=0.5*gPart[noffset[i]+j].GetMass()*v2;
#ifdef GASON
            Ti+=gPart[noffset[i]+j].GetMass()*gPart[noffset[i]+j].GetU();
#endif
#endif
            totT+=Ti;
            gPart[noffset[i]+j].SetDensity(opt.uinfo.Eratio*Ti+gPart[noffset[i]+j].GetPotential());
            Efrac+=(Ti+gPart[noffset[i]+j].GetPotential()<0);
            if(maxE<gPart[noffset[i]+j].GetDensity()) maxE=gPart[noffset[i]+j].GetDensity();
        }
        Efrac/=(Double_t)numingroup[i];
        //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
        pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+2);
        Eplusflag=new int[numingroup[i]];
        nEplusid=new Int_t[pqsize];
        if (opt.uinfo.unbindtype==USYSANDPART) {
            if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
            else unbindcheck=false;
        }
        else if (opt.uinfo.unbindtype==UPART) {
            if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
            else unbindcheck=false;
        }
        if (unbindcheck) {
            for (j=0;j<numingroup[i];j++)Eplusflag[j]=0;
            pq=new PriorityQueue(pqsize);
            for (j=0;j<pqsize;j++) pq->Push(j,gPart[noffset[i]+j].GetDensity());
            for (j=pqsize;j<numingroup[i];j++) if (gPart[noffset[i]+j].GetDensity()>gPart[noffset[i]+pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[noffset[i]+j].GetDensity());}
            nEplus=0;
            //if just looking at particle then add to removal list till energy >0
            if (opt.uinfo.unbindtype==UPART) {
                for (j=0;j<pqsize;j++) {
                    if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
            else if (opt.uinfo.unbindtype==USYSANDPART) {
                nEfrac=0;
                if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                for (j=0;j<pqsize;j++) {
                    if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                    else break;
                }
            }
            delete pq;
        }
        while(unbindcheck)
        {
            iunbindflag++;
            //first correct for removal of all least bound particle
            double temp=1.0/gmass[i], temp2=0.;
            if (opt.uinfo.cmvelreftype==CMVELREF) {
                for (j=0;j<nEplus;j++) {
                    for (k=0;k<3;k++) cmvel[i][k]-=gPart[noffset[i]+nEplusid[j]].GetVelocity(k)*gPart[noffset[i]+nEplusid[j]].GetMass()*temp;
                    temp2+=gPart[noffset[i]+nEplusid[j]].GetMass();
                }
                temp=gmass[i]/(gmass[i]-temp2);
                for (k=0;k<3;k++) cmvel[i][k]*=temp;
                gmass[i]-=temp2;
            }
            else {
                for (j=0;j<nEplus;j++) gmass[i]-=gPart[noffset[i]+nEplusid[j]].GetMass();
            }
            //if ignore the background then adjust the potential energy of the particles
            if (opt.uinfo.bgpot==0) {
                for (k=0;k<nEplus;k++) {
                    totV[i]-=0.5*gPart[noffset[i]+nEplusid[k]].GetPotential();
                    for (j=0;j<numingroup[i];j++) {
                        if (j!=nEplusid[k]) {
                            r2=0.;for (n=0;n<3;n++) r2+=pow(gPart[noffset[i]+nEplusid[k]].GetPosition(n)-gPart[noffset[i]+j].GetPosition(n),2.0);
                            r2+=eps2;
                            r2=1.0/sqrt(r2);
#ifdef NOMASS
                            poti=gPart[noffset[i]+j].GetPotential()+opt.G*(gPart[noffset[i]+nEplusid[k]].GetMass()*gPart[noffset[i]+j].GetMass())*r2*mv2;
#else
                            poti=gPart[noffset[i]+j].GetPotential()+opt.G*(gPart[noffset[i]+nEplusid[k]].GetMass()*gPart[noffset[i]+j].GetMass())*r2;
#endif
                            gPart[noffset[i]+j].SetPotential(poti);
                        }
                    }
                }
            }
            //remove particles with positive energy
            for (j=0;j<nEplus;j++) pfof[gPart[noffset[i]+nEplusid[j]].GetPID()]=0;
            k=numingroup[i]-1;
            for (j=0;j<nEplus;j++) if (nEplusid[j]<numingroup[i]-nEplus) {
                while(Eplusflag[k]==1)k--;
                Ptemp=gPart[noffset[i]+nEplusid[j]];
                gPart[noffset[i]+nEplusid[j]]=gPart[noffset[i]+k];
                gPart[noffset[i]+k]=Ptemp;
                Eplusflag[nEplusid[j]]=0;
                k--;
            }
            numingroup[i]-=nEplus;

            if (nEplus>=0.1*pqsize+0.5) {

            //recalculate kinetic energies since cmvel has changed
            totT=0;
            maxE=-MAXVALUE;
            for (j=0;j<numingroup[i];j++) {
                v2=0.0;for (k=0;k<3;k++) v2+=pow(gPart[noffset[i]+j].GetVelocity(k)-cmvel[i][k],2.0);
#ifdef NOMASS
                Ti=0.5*gPart[noffset[i]+j].GetMass()*v2*opt.MassValue;
#ifdef GASON
                Ti+=opt.MassValue*gPart[noffset[i]+j].GetU();
#endif
#else
                Ti=0.5*gPart[noffset[i]+j].GetMass()*v2;
#ifdef GASON
                Ti+=gPart[noffset[i]+j].GetMass()*gPart[noffset[i]+j].GetU();
#endif
#endif
                totT+=Ti;
                gPart[noffset[i]+j].SetDensity(opt.uinfo.Eratio*Ti+gPart[noffset[i]+j].GetPotential());
                Efrac+=(Ti+gPart[noffset[i]+j].GetPotential()<0);
                    if(maxE<gPart[noffset[i]+j].GetDensity()) maxE=gPart[noffset[i]+j].GetDensity();
            }
            Efrac/=(Double_t)numingroup[i];
            if (opt.uinfo.unbindtype==USYSANDPART) {
                if(((Efrac<opt.uinfo.minEfrac)||(maxE>0))&&(numingroup[i]>=opt.MinSize)) unbindcheck=true;
                else unbindcheck=false;
            }
            else if (opt.uinfo.unbindtype==UPART) {
                if ((maxE>0)&&(numingroup[i]>=opt.MinSize))unbindcheck=true;
                else unbindcheck=false;
            }
            if (unbindcheck) {
                //determine if any particle  number of particle with positive energy upto opt.uinfo.maxunbindfrac*numingroup+1
                pqsize=(Int_t)(opt.uinfo.maxunbindfrac*numingroup[i]+1);
                pq=new PriorityQueue(pqsize);
                for (j=0;j<pqsize;j++) pq->Push(j,gPart[noffset[i]+j].GetDensity());
                for (j=pqsize;j<numingroup[i];j++) if (gPart[noffset[i]+j].GetDensity()>gPart[noffset[i]+pq->TopQueue()].GetDensity()) {pq->Pop();pq->Push(j,gPart[noffset[i]+j].GetDensity());}
                nEplus=0;
                //if just looking at particle then add to removal list till energy >0
                if (opt.uinfo.unbindtype==UPART) {
                    for (j=0;j<pqsize;j++) {
                        if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                //otherwise, remove all positive energies and also if Efrac< minEfrac, keep adding to removal list
                else if (opt.uinfo.unbindtype==USYSANDPART) {
                    nEfrac=0;
                    if (Efrac<opt.uinfo.minEfrac) nEfrac=(opt.uinfo.minEfrac-Efrac)*numingroup[i];
                    for (j=0;j<pqsize;j++) {
                        if (gPart[noffset[i]+pq->TopQueue()].GetDensity()>0 || nEplus<nEfrac) {nEplusid[nEplus++]=pq->TopQueue();Eplusflag[pq->TopQueue()]=1;pq->Pop();}
                        else break;
                    }
                }
                delete pq;
            }
            }
            else unbindcheck=false;
        }
        //if group too small remove entirely
        if (numingroup[i]<opt.MinSize) {
            iunbindflag++;
            for (j=0;j<numingroup[i];j++) pfof[gPart[noffset[i]+j].GetPID()]=0;
            numingroup[i]=0;
        }
        delete[] nEplusid;
        delete[] Eplusflag;
    }
#ifdef USEOPENMP
}
#endif
    for (i=1;i<=numgroups;i++) if (numingroup[i]==0) ng--;
    delete[] cmvel;
    numgroups=ng;
    //return if any unbinding done indicating groups have been reordered
    if (iunbindflag) return 1;
    else return 0;
}

/// Calculates the gravitational potential using a kd-tree and monopole expansion
///\todo need ewald correction for periodic systems and also use more than monopole.
void Potential(Options &opt, Int_t nbodies, Particle *Part, Double_t *potV)
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

    cout<<"Calculating potentials ..."<<endl;
    //otherwise use tree tree gravity calculation
    //here openmp is per group since each group is large
    //to make this memory efficient really need just KDTree that uses Coordinates
    tree=new KDTree(Part,nbodies,opt.uinfo.BucketSize,tree->TPHYS,tree->KEPAN,1000,1);
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
        potV[Part[j].GetID()]=0.;
        ntreecell=nleafcell=0;
        Coordinate xpos(Part[j].GetPosition());
        MarkCell(npomp[tid],marktreecell[tid], markleafcell[tid],ntreecell,nleafcell,r2val[tid],opt.uinfo.BucketSize, cR2max, cellcm, cmtot, xpos, eps2);
        for (k=0;k<ntreecell;k++) {
          potV[Part[j].GetID()]+=-Part[j].GetMass()*r2val[tid][k];
        }
        for (k=0;k<nleafcell;k++) {
            for (l=start[markleafcell[tid][k]];l<end[markleafcell[tid][k]];l++) {
                if (j!=l) {
                    r2=0.;for (n=0;n<3;n++) r2+=pow(Part[j].GetPosition(n)-Part[l].GetPosition(n),(Double_t)2.0);
                    r2+=eps2;
                    r2=1.0/sqrt(r2);
                    potV[Part[j].GetID()]+=-(Part[j].GetMass()*Part[l].GetMass())*r2;
                }
            }
        }
        potV[Part[j].GetID()]*=opt.G;
    }
#ifdef USEOPENMP
}
#endif
    delete tree;
    cout<<"Done\n";
}

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
    delete[] start;
    delete[] end;
    delete[] cmtot;
    delete[] cBmax;
    delete[] cR2max;
    delete[] cellcm;
    delete[] nodelist;
    for (j=0;j<nthreads;j++) {delete[] marktreecell[j]; delete[] markleafcell[j]; delete[] r2val[j];}
    delete[] marktreecell;
    delete[] markleafcell;
    delete[] r2val;
    delete[] npomp;
}
