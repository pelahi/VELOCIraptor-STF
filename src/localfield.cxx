/*! \file localfield.cxx
 *  \brief this file contains routines to calculate the local velocity density function for each particle
 */

//--  Local Velocity density routines

#include "stf.h"
#include "swiftinterface.h"

/*! Calculates the local velocity density function for each particle using a kernel technique
    There are two approaches to getting this local quantity \n
    1) From a large set of nearest physical neighbours use a smaller subset of nearest velocity neighbours \n
    2) Or calculate phase space density and physical density and divide phase-space density by physical density. This is far more computationally expensive and may not enhance velocity clustering
    and is just present for testing purposes.
    \todo velocity density function is NOT mass weighted. Might want to alter this.
    \todo there is a seg fault memory error when searching for NN in large sims using \em SINGLEPRECISION flag. I don't know why.
*/
void GetVelocityDensity(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree)
{
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    auto time1 = MyGetTime();
    cout<<ThisTask<<": Get local velocity density"<<endl;
    if (opt.iverbose) {
        cout<<ThisTask<<" "<<"Using the following parameters to calculate velocity density using sph kernel: ";
        cout<<ThisTask<<" "<<"(Nse,Nv)="<<opt.Nsearch<<","<<opt.Nvel<<endl;
        cout<<ThisTask<<" "<<"Get velocity density using a subset of nearby physical or phase-space neighbours"<<endl;
        if (tree == NULL) cout<<ThisTask<<" Building Tree first in (x) space to get local velocity density"<<endl;
    }
#ifdef HALOONLYDEN
    GetVelocityDensityHaloOnlyDen(opt, nbodies, Part, tree);
#else
    if (opt.iLocalVelDenApproxCalcFlag>0) GetVelocityDensityApproximative(opt, nbodies, Part, tree);
    else GetVelocityDensityExact(opt, nbodies, Part, tree);
#endif
    cout<<ThisTask<<": finished calculation in "<<MyElapsedTime(time1)<<endl;
}

void GetVelocityDensityOld(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree)
{
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
Int_t i,j,k;
    int nthreads;
    int tid,id,pid,pid2,itreeflag=0;
    Double_t v2;
    Double_t *period=NULL;
    ///\todo alter period so arbitrary dimensions
    if (opt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=opt.p;
    }
    //if using mpi run NN search store largest distance for each particle so that export list can be built.
    //if calculating using only particles IN a structure,
#ifndef HALOONLYDEN
#ifdef USEMPI
    Int_t nimport;
    Double_t *maxrdist=new Double_t[nbodies];
    Double_t *weight;
    Int_t *nnids,*nnidsneighbours;
    Double_t *nnr2, *nnr2neighbours;
    PriorityQueue *pqx, *pqv;
#ifndef USEOPENMP
    nthreads=1;
#else
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

    auto time2 = MyGetTime();
    //In loop determine if particles NN search radius overlaps another mpi threads domain.
    //If not, then proceed as usually to determine velocity density.
    //If so, do not calculate local velocity density and set its velocity density to -1 as a flag
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k,tid,id,v2,nnids,nnr2,nnidsneighbours,nnr2neighbours,weight,pqx,pqv)
{
#endif
    nnids=new Int_t[opt.Nsearch];
    nnr2=new Double_t[opt.Nsearch];
    weight=new Double_t[opt.Nvel];
    pqx=new PriorityQueue(opt.Nsearch);
    pqv=new PriorityQueue(opt.Nvel);
#ifdef USEOPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i=0;i<nbodies;i++) {
        //if strucden compile flag set then only calculate velocity density for particles in groups
#ifdef STRUCDEN
        if (Part[i].GetType()>0) {
        //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
        if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
        //otherwise distinction must be made so that only base calculation on dark matter particles
        else tree->FindNearestCriterion(i,FOFPositivetypes,NULL,nnids,nnr2,opt.Nsearch);
#else
        tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
#endif
        //once NN set is found, store maxrdist and see if particle's search radius overlaps with another mpi domain
        maxrdist[i]=sqrt(nnr2[opt.Nsearch-1]);
        bool ioverlap;
        if (opt.impiusemesh) {
            ioverlap = (MPISearchForOverlapUsingMesh(opt,Part[i],maxrdist[i])==0);
        }
        else
        {
            ioverlap = (MPISearchForOverlap(Part[i],maxrdist[i])==0);
        }
        if (ioverlap) {
            for (j=0;j<opt.Nvel;j++) {
                pqv->Push(-1, MAXVALUE);
                weight[j]=1.0;
            }
            for (j=0;j<opt.Nsearch;j++) {
                v2=0;
                id=nnids[j];
                for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-Part[id].GetVelocity(k))*(Part[i].GetVelocity(k)-Part[id].GetVelocity(k));
                if (v2 < pqv->TopPriority()){
                    pqv->Pop();
                    pqv->Push(id, v2);
                }
            }
            Part[i].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv, weight));
        }
        else Part[i].SetDensity(-1.0);
#ifdef STRUCDEN
        }
        else maxrdist[i]=0.0;
#endif
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] weight;
    delete pqx;
    delete pqv;
#ifdef USEOPENMP
}
#endif
    if (opt.iverbose) cout<<ThisTask<<" finished local calculation in "<<MyElapsedTime(time2)<<endl;
    time2=MyGetTime();

    //determines export AND import numbers

    if (opt.impiusemesh) MPIGetNNExportNumUsingMesh(opt, nbodies, Part, maxrdist);
    else MPIGetNNExportNum(nbodies, Part, maxrdist);
    NNDataIn = new nndata_in[NExport];
    NNDataGet = new nndata_in[NImport];
    //build the exported particle list using NNData structures
    if (opt.impiusemesh) MPIBuildParticleNNExportListUsingMesh(opt, nbodies, Part, maxrdist);
    else MPIBuildParticleNNExportList(nbodies, Part, maxrdist);
    MPIGetNNImportNum(nbodies, tree, Part, (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)));
    PartDataIn = new Particle[NExport];
    PartDataGet = new Particle[NImport];
    MPI_Barrier(MPI_COMM_WORLD);
    //run search on exported particles and determine which local particles need to be exported back (or imported)
    nimport=MPIBuildParticleNNImportList(opt, nbodies, tree, Part, (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)));
    int nimportsearch=opt.Nsearch;
    if (nimportsearch>nimport) nimportsearch=nimport;
    if (opt.iverbose) cout<<ThisTask<<" Searching particles in other domains "<<nimport<<endl;
    //now with imported particle list and local particle list can run proper NN search
    //first build neighbouring tree
    KDTree *treeneighbours=NULL;
    if (nimport>0) treeneighbours=new KDTree(PartDataGet,nimport,1,tree->TPHYS,tree->KEPAN,100,0,0,0,period);
    //then run search
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k,tid,pid,pid2,v2,nnids,nnr2,nnidsneighbours,nnr2neighbours,weight,pqx,pqv)
{
#endif
    nnids=new Int_t[opt.Nsearch];
    nnr2=new Double_t[opt.Nsearch];
    nnidsneighbours=new Int_t[nimportsearch];
    nnr2neighbours=new Double_t[nimportsearch];
    weight=new Double_t[opt.Nvel];
    pqx=new PriorityQueue(opt.Nsearch);
    pqv=new PriorityQueue(opt.Nvel);
#ifdef USEOPENMP
#pragma omp for
#endif
    for (i=0;i<nbodies;i++) {
#ifdef STRUCDEN
        if (Part[i].GetType()>0) {
#endif
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        pid=Part[i].GetID();
        if (Part[i].GetDensity()==-1) {
            //search trees

            //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
            if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
            //otherwise distinction must be made so that only base calculation on dark matter particles
            else tree->FindNearestCriterion(i,FOFPositivetypes,NULL,nnids,nnr2,opt.Nsearch);

            //fill priority queue, where queue value for neighbours is offset by local particle number
            for (j = 0; j <opt.Nsearch; j++) pqx->Push(-1, MAXVALUE);
            for (j=0;j<opt.Nsearch;j++) {
                if (nnr2[j] < pqx->TopPriority()){
                    pqx->Pop();
                    pqx->Push(nnids[j], nnr2[j]);
                }
            }
            //now search the export particle list and fill appropriately
            if (nimport>0) {
                Coordinate x(Part[i].GetPosition());
                treeneighbours->FindNearestPos(x,nnidsneighbours,nnr2neighbours,nimportsearch);
                for (j=0;j<nimportsearch;j++) {
                    if (nnr2neighbours[j] < pqx->TopPriority()){
                        pqx->Pop();
                        pqx->Push(nnidsneighbours[j]+nbodies, nnr2neighbours[j]);
                    }
                }
            }
            for (j=0;j<opt.Nvel;j++) {
                pqv->Push(-1, MAXVALUE);
                weight[j]=1.0;
            }
            for (j=0;j<opt.Nsearch;j++) {
                v2=0;
                if (pqx->TopQueue()<nbodies) {
                    pid2=pqx->TopQueue();
                    for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-Part[pid2].GetVelocity(k))*(Part[i].GetVelocity(k)-Part[pid2].GetVelocity(k));
                }
                else {
                    pid2=pqx->TopQueue()-nbodies;
                    for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-PartDataGet[pid2].GetVelocity(k))*(Part[i].GetVelocity(k)-PartDataGet[pid2].GetVelocity(k));
                }
                if (v2 < pqv->TopPriority()){
                    pqv->Pop();
                    pqv->Push(pqx->TopQueue(), v2);
                }
                pqx->Pop();
            }
            //and now calculate velocity density function
            Part[i].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv, weight));
        }
#ifdef STRUCDEN
        }
#endif
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] nnidsneighbours;
    delete[] nnr2neighbours;
    delete[] weight;
    delete pqx;
    delete pqv;
#ifdef USEOPENMP
}
#endif
    //free memory
    if (itreeflag) delete tree;
    if (nimport>0) delete treeneighbours;
    delete[] PartDataIn;
    delete[] PartDataGet;
    delete[] NNDataIn;
    delete[] NNDataGet;
    if(opt.iverbose) cout<<ThisTask<<" finished other domain search "<<MyElapsedTime(time2)<<endl;
#else
    //NO MPI invoked
#ifndef USEOPENMP
    for (i=0;i<nbodies;i++) {
#ifdef STRUCDEN
        if (Part[i].GetType()>0)
        {
#endif
            Part[i].SetDensity(tree->CalcVelDensityParticle(i,opt.Nvel,opt.Nsearch));
#ifdef STRUCDEN
        }
#endif
    }
#else
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
    Int_t *fracdone=new Int_t[nthreads];
    Int_t *fraclim=new Int_t[nthreads];
    int minamount=(double)nbodies/(double)nthreads*0.01+1;
    Int_t addamount=(double)nbodies/(double)nthreads*0.1+1;
    for (int j=0;j<nthreads;j++){fracdone[j]=0;fraclim[j]=addamount;}
    Int_t *nnids;
    Double_t *nnr2;
    Double_t *weight;
    PriorityQueue **pqx, **pqv;

    nnids=new Int_t[nthreads*opt.Nsearch];
    nnr2=new Double_t[nthreads*opt.Nsearch];
    pqx=new PriorityQueue*[nthreads];
    pqv=new PriorityQueue*[nthreads];
    weight=new Double_t[nthreads*opt.Nvel];
    for (j=0;j<nthreads;j++) {
        //nnids[j]=new Int_t[opt.Nsearch];
        //nnr2[j]=new Double_t[opt.Nsearch];
        pqx[j]=new PriorityQueue(opt.Nsearch);
        pqv[j]=new PriorityQueue(opt.Nvel);
    }
#pragma omp parallel default(shared) \
private(i,tid)
{
#pragma omp for schedule(dynamic,minamount) nowait
    for (i=0;i<nbodies;i++) {
#ifdef STRUCDEN
        if (Part[i].GetType()>0) {
#endif
        tid=omp_get_thread_num();
        //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
        if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) Part[i].SetDensity(tree->CalcVelDensityParticle(i,opt.Nvel,opt.Nsearch,1,pqx[tid],pqv[tid],&nnids[tid*opt.Nsearch],&nnr2[tid*opt.Nsearch]));
        //otherwise distinction must be made so that only base calculation on dark matter particles
        else {
            tree->FindNearestCriterion(i,FOFPositivetypes,NULL,&nnids[tid*opt.Nsearch],&nnr2[tid*opt.Nsearch],opt.Nsearch);
            for (j=0;j<opt.Nvel;j++) {
                pqv[tid]->Push(-1, MAXVALUE);
                weight[j+tid*opt.Nvel]=1.0;
            }
            for (j=0;j<opt.Nsearch;j++) {
                v2=0;
                id=nnids[j+tid*opt.Nsearch];
                for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-Part[id].GetVelocity(k))*(Part[i].GetVelocity(k)-Part[id].GetVelocity(k));
                if (v2 < pqv[tid]->TopPriority()){
                    pqv[tid]->Pop();
                    pqv[tid]->Push(id, v2);
                }
            }
            Part[i].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv[tid], &weight[tid*opt.Nvel]));
        }

        fracdone[tid]++;
        if (opt.iverbose) if (fracdone[tid]>fraclim[tid]) {printf("Task %d done %e of its share \n",tid,fracdone[tid]/((double)nbodies/(double)nthreads));fraclim[tid]+=addamount;}
#ifdef STRUCDEN
        }
#endif
    }
}
#endif
    if (itreeflag) delete tree;
#ifdef USEOPENMP
    for (j=0;j<nthreads;j++) {
        delete pqx[j];
        delete pqv[j];
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] pqx;
    delete[] pqv;
    delete[] fracdone;
    delete[] fraclim;
#endif
#endif
#else

    //start halo only density calculations, where particles are localized to single mpi domain
    nthreads=1;
#ifdef USEOPENMP
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

    Int_t *fracdone=new Int_t[nthreads];
    Int_t *fraclim=new Int_t[nthreads];
    int minamount=(double)nbodies/(double)nthreads*0.01+1;
    Int_t addamount=(double)nbodies/(double)nthreads*0.1+1;
    for (j=0;j<nthreads;j++){fracdone[j]=0;fraclim[j]=addamount;}
    Int_t *nnids;
    Double_t *nnr2;
    PriorityQueue **pqx, **pqv;

    nnids=new Int_t[nthreads*opt.Nsearch];
    nnr2=new Double_t[nthreads*opt.Nsearch];
    pqx=new PriorityQueue*[nthreads];
    pqv=new PriorityQueue*[nthreads];
    for (j=0;j<nthreads;j++) {
        //nnids[j]=new Int_t[opt.Nsearch];
        //nnr2[j]=new Double_t[opt.Nsearch];
        pqx[j]=new PriorityQueue(opt.Nsearch);
        pqv[j]=new PriorityQueue(opt.Nvel);
    }
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid)
{
#pragma omp for schedule(dynamic,minamount) nowait
#endif
    for (i=0;i<nbodies;i++) {
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        Part[i].SetDensity(tree->CalcVelDensityParticle(i,opt.Nvel,opt.Nsearch,1,pqx[tid],pqv[tid],&nnids[tid*opt.Nsearch],&nnr2[tid*opt.Nsearch]));
        fracdone[tid]++;
        if (opt.iverbose) if (fracdone[tid]>fraclim[tid]) {printf("Task %d done %e of its share \n",tid,fracdone[tid]/((double)nbodies/(double)nthreads));fraclim[tid]+=addamount;}
    }
#ifdef USEOPENMP
}
#endif
    for (j=0;j<nthreads;j++) {
        //delete[] nnids[j];
        //delete[] nnr2[j];
        delete pqx[j];
        delete pqv[j];
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] pqx;
    delete[] pqv;
    delete[] fracdone;
    delete[] fraclim;
    if (itreeflag) delete tree;
#endif
    if (period!=NULL) delete[] period;
}

//start halo only density calculations, where particles are localized to single halo
void GetVelocityDensityHaloOnlyDen(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree)
{
    Int_t i,j;
    int nthreads;
    int tid;
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
#endif
    nthreads=1;
#ifdef USEOPENMP
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

    tree=new KDTree(Part,nbodies,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,0,0,0);
    Int_t *nnids;
    Double_t *nnr2;
    PriorityQueue **pqx, **pqv;

    nnids=new Int_t[nthreads*opt.Nsearch];
    nnr2=new Double_t[nthreads*opt.Nsearch];
    pqx=new PriorityQueue*[nthreads];
    pqv=new PriorityQueue*[nthreads];
    for (j=0;j<nthreads;j++) {
        pqx[j]=new PriorityQueue(opt.Nsearch);
        pqv[j]=new PriorityQueue(opt.Nvel);
    }

    //get memory useage
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=1));

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tid)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<nbodies;i++) {
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        Part[i].SetDensity(tree->CalcVelDensityParticle(i,opt.Nvel,opt.Nsearch,1,pqx[tid],pqv[tid],&nnids[tid*opt.Nsearch],&nnr2[tid*opt.Nsearch]));
    }
#ifdef USEOPENMP
}
#endif
    for (j=0;j<nthreads;j++) {
        delete pqx[j];
        delete pqv[j];
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] pqx;
    delete[] pqv;
    delete tree;
}

///Exact calculation of velocity density at a particle's position
void GetVelocityDensityExact(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree)
{
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    if (opt.iverbose) cout<<ThisTask<<" Calculating the local velocity density by finding EXACT nearest physical neighbours to particles"<<endl;
    Int_t i,j,k;
    int nthreads;
    int tid,id,pid,pid2,itreeflag=0;
    Double_t v2;
    Double_t *period=NULL;
    ///\todo alter period so arbitrary dimensions
    if (opt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=opt.p;
    }
    //if using mpi run NN search store largest distance for each particle so that export list can be built.
    //if calculating using only particles IN a structure,
    Double_t *weight;
    Int_t *nnids;
    Double_t *nnr2;
    PriorityQueue *pqx, *pqv;

#ifdef USEMPI
    Int_t nimport;
    Int_t *nnidsneighbours;
    Double_t *nnr2neighbours;
    Double_t *maxrdist=NULL;
    maxrdist = new Double_t[nbodies];
    for (i=0;i<nbodies;i++) maxrdist[i]=0;
#endif

#ifndef USEOPENMP
    nthreads=1;
#else
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

    auto time2=MyGetTime();
    //get memory useage
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=1));

    //In loop determine if particles NN search radius overlaps another mpi threads domain.
    //If not, then proceed as usually to determine velocity density.
    //If so, do not calculate local velocity density and set its velocity density to -1 as a flag
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k,tid,id,v2,nnids,nnr2,weight,pqv)
{
#endif
    nnids=new Int_t[opt.Nsearch];
    nnr2=new Double_t[opt.Nsearch];
    weight=new Double_t[opt.Nvel];
    pqv=new PriorityQueue(opt.Nvel);
#ifdef USEOPENMP
#pragma omp for schedule(dynamic)
#endif
    for (i=0;i<nbodies;i++) {
        //if strucden compile flag set then only calculate velocity density for particles in groups
#ifdef STRUCDEN
        if (Part[i].GetType()<=0) continue;
        //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
        if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
        //otherwise distinction must be made so that only base calculation on dark matter particles
        else tree->FindNearestCriterion(i,FOFPositivetypes,NULL,nnids,nnr2,opt.Nsearch);
#else
        tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
#endif
#ifdef USEMPI
        if (opt.iLocalVelDenApproxCalcFlag==0 && NProcs>1) {
            //once NN set is found, store maxrdist and see if particle's search radius overlaps with another mpi domain
            maxrdist[i]=sqrt(nnr2[opt.Nsearch-1]);
            bool ioverlap;

            if (opt.impiusemesh) ioverlap = (MPISearchForOverlapUsingMesh(opt,Part[i],maxrdist[i])!=0);
            else ioverlap = (MPISearchForOverlap(Part[i],maxrdist[i])!=0);
            if (ioverlap) {
                Part[i].SetDensity(-1.0);
                continue;
            }
            maxrdist[i]=0.0;
        }
#endif
        for (j=0;j<opt.Nvel;j++) {
            pqv->Push(-1, MAXVALUE);
            weight[j]=1.0;
        }
        for (j=0;j<opt.Nsearch;j++) {
            v2=0;
            id=nnids[j];
            for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-Part[id].GetVelocity(k))*(Part[i].GetVelocity(k)-Part[id].GetVelocity(k));
            if (v2 < pqv->TopPriority()){
                pqv->Pop();
                pqv->Push(id, v2);
            }
        }
        Part[i].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv, weight));
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] weight;
    delete pqv;
#ifdef USEOPENMP
}
#endif

#ifdef USEMPI
    if (NProcs >1 && opt.iLocalVelDenApproxCalcFlag==0) {
    if (opt.iverbose) cout<<ThisTask<<" finished local calculation in "<<MyElapsedTime(time2)<<endl;
    time2=MyGetTime();
    //determines export AND import numbers
    if (opt.impiusemesh) MPIGetNNExportNumUsingMesh(opt, nbodies, Part, maxrdist);
    else MPIGetNNExportNum(nbodies, Part, maxrdist);
    NNDataIn = new nndata_in[NExport];
    NNDataGet = new nndata_in[NImport];
    //build the exported particle list using NNData structures
    if (opt.impiusemesh) MPIBuildParticleNNExportListUsingMesh(opt, nbodies, Part, maxrdist);
    else MPIBuildParticleNNExportList(nbodies, Part, maxrdist);
    MPIGetNNImportNum(nbodies, tree, Part, (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)));
    PartDataIn = new Particle[NExport];
    PartDataGet = new Particle[NImport];
    MPI_Barrier(MPI_COMM_WORLD);
    //run search on exported particles and determine which local particles need to be exported back (or imported)
    nimport=MPIBuildParticleNNImportList(opt, nbodies, tree, Part,(!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)));
    int nimportsearch=opt.Nsearch+1;
    if (nimportsearch>nimport) nimportsearch=nimport;
    if (opt.iverbose) cout<<ThisTask<<" Searching particles in other domains "<<nimport<<endl;
    //now with imported particle list and local particle list can run proper NN search
    //first build neighbouring tree
    KDTree *treeneighbours=NULL;
    if (nimport>0) treeneighbours=new KDTree(PartDataGet,nimport,1,tree->TPHYS,tree->KEPAN,100,0,0,0,period);

    //get memory useage
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=1));

    //then run search
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k,tid,pid,pid2,v2,nnids,nnr2,nnidsneighbours,nnr2neighbours,weight,pqx,pqv)
{
#endif
    nnids=new Int_t[opt.Nsearch];
    nnr2=new Double_t[opt.Nsearch];
    nnidsneighbours=new Int_t[nimportsearch];
    nnr2neighbours=new Double_t[nimportsearch];
    weight=new Double_t[opt.Nvel];
    pqx=new PriorityQueue(opt.Nsearch);
    pqv=new PriorityQueue(opt.Nvel);
#ifdef USEOPENMP
#pragma omp for
#endif
    for (i=0;i<nbodies;i++) {
#ifdef STRUCDEN
        if (Part[i].GetType()<=0) continue;
#endif

#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        pid=Part[i].GetID();
        if (Part[i].GetDensity()==-1) {
            //search trees

            //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
            if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
            //otherwise distinction must be made so that only base calculation on dark matter particles
            else tree->FindNearestCriterion(i,FOFPositivetypes,NULL,nnids,nnr2,opt.Nsearch);

            //fill priority queue, where queue value for neighbours is offset by local particle number
            for (j = 0; j <opt.Nsearch; j++) pqx->Push(-1, MAXVALUE);
            for (j=0;j<opt.Nsearch;j++) {
                if (nnr2[j] < pqx->TopPriority()){
                    pqx->Pop();
                    pqx->Push(nnids[j], nnr2[j]);
                }
            }
            //now search the export particle list and fill appropriately
            if (nimport>0) {
                Coordinate x(Part[i].GetPosition());
                treeneighbours->FindNearestPos(x,nnidsneighbours,nnr2neighbours,nimportsearch);
                for (j=0;j<nimportsearch;j++) {
                    if (nnr2neighbours[j] < pqx->TopPriority()){
                        pqx->Pop();
                        pqx->Push(nnidsneighbours[j]+nbodies, nnr2neighbours[j]);
                    }
                }
            }
            for (j=0;j<opt.Nvel;j++) {
                pqv->Push(-1, MAXVALUE);
                weight[j]=1.0;
            }
            for (j=0;j<opt.Nsearch;j++) {
                v2=0;
                if (pqx->TopQueue()<nbodies) {
                    pid2=pqx->TopQueue();
                    for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-Part[pid2].GetVelocity(k))*(Part[i].GetVelocity(k)-Part[pid2].GetVelocity(k));
                }
                else {
                    pid2=pqx->TopQueue()-nbodies;
                    for (k=0;k<3;k++) v2+=(Part[i].GetVelocity(k)-PartDataGet[pid2].GetVelocity(k))*(Part[i].GetVelocity(k)-PartDataGet[pid2].GetVelocity(k));
                }
                if (v2 < pqv->TopPriority()){
                    pqv->Pop();
                    pqv->Push(pqx->TopQueue(), v2);
                }
                pqx->Pop();
            }
            //and now calculate velocity density function
            Part[i].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv, weight));
        }
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] nnidsneighbours;
    delete[] nnr2neighbours;
    delete[] weight;
    delete pqx;
    delete pqv;
#ifdef USEOPENMP
}
#endif
    //free memory
    if (nimport>0) delete treeneighbours;
    delete[] PartDataIn;
    delete[] PartDataGet;
    delete[] NNDataIn;
    delete[] NNDataGet;
    if(opt.iverbose) cout<<ThisTask<<" finished other domain search "<<MyElapsedTime(time2)<<endl;
    }
#endif
    if (itreeflag) delete tree;
    if (period!=NULL) delete[] period;
}

void GetVelocityDensityApproximative(Options &opt, const Int_t nbodies, Particle *Part, KDTree *tree)
{
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
#endif
    if (opt.iverbose) cout<<ThisTask<<" Calculating the local velocity density by finding APPROXIMATIVE nearest physical neighbour search for each particle "<<endl;
    int nthreads;
    int id,pid2,itreeflag=0;
    Double_t v2;
    Int_t nprocessed=0, ntot=0;
    ///\todo alter period so arbitrary dimensions
    Double_t *period=NULL;
    if (opt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=opt.p;
    }
    //if using mpi run NN search store largest distance for each particle so that export list can be built.
    //if calculating using only particles IN a structure,
    Int_t nimport;
    Double_t *maxrdist=NULL;
    Double_t *weight;
    Int_t *nnids,*nnidsneighbours;
    Double_t *nnr2, *nnr2neighbours;
    PriorityQueue *pqx, *pqv;
    Particle *Pval;
#ifndef USEOPENMP
    nthreads=1;
#else
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif

    auto time2=MyGetTime();
    //only build tree if necessary
    if (tree==NULL) {
        itreeflag=1;
        tree=new KDTree(Part,nbodies,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,0,0,0,period);
    }
    //In loop determine if particles NN search radius overlaps another mpi threads domain.
    //If not, then proceed as usually to determine velocity density.
    //If so, do not calculate local velocity density and set its velocity density to -1 as a flag
    //idea is to use approximative near neighbour search using the centre-of-mass of the bucket
    //and the size of the bucket to get particle list.
    //if MPI is used, then also use size of bucket to determine rough maximum search distance
    //and whether other mpi domains need to be searched.

    //first get all local leaf nodes;
    Int_t numleafnodes = tree->GetNumLeafNodes();
    Node *node;
    vector<leaf_node_info> leafnodes(numleafnodes);
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

    //get memory useage
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=1));

#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
#pragma omp for schedule(dynamic)
#endif
    for (auto i=0;i<numleafnodes;i++)
    {
        leafnodes[i].num=0;
        leafnodes[i].cm[0]=leafnodes[i].cm[1]=leafnodes[i].cm[2]=0;
#ifdef USEMPI
        leafnodes[i].searchdist = 0;
#endif
        for (auto j=leafnodes[i].istart;j<leafnodes[i].iend;j++)
        {
#ifdef STRUCDEN
            if (Part[j].GetType()<=0) continue;
#endif
            leafnodes[i].num++;
            for (auto k=0;k<3;k++) leafnodes[i].cm[k] += Part[j].GetPosition(k);
        }
        if (leafnodes[i].num == 0) continue;
        for (auto k=0;k<3;k++) leafnodes[i].cm[k]/= (float)leafnodes[i].num;
        for (auto j=leafnodes[i].istart;j<leafnodes[i].iend;j++)
        {
            double r2 = 0;
            for (auto k=0;k<3;k++) r2 += pow(Part[j].GetPosition(k)-leafnodes[i].cm[k],2.0);
            if (r2>leafnodes[i].size) leafnodes[i].size=r2;
        }
        leafnodes[i].size=sqrt(leafnodes[i].size);
#ifdef USEMPI
        leafnodes[i].searchdist=leafnodes[i].size*pow(opt.Nsearch/(float)leafnodes[i].numtot,1.0/3.0)*1.2;
#endif
    }
#ifdef USEOPENMP
}
#endif

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(id,v2,nnids,nnr2,weight,pqv)
{
#endif
    nnids=new Int_t[opt.Nsearch];
    nnr2=new Double_t[opt.Nsearch];
    weight=new Double_t[opt.Nvel];
    pqv=new PriorityQueue(opt.Nvel);
#ifdef USEOPENMP
#pragma omp for schedule(dynamic) \
reduction(+:nprocessed,ntot)
#endif
    for (auto i=0;i<numleafnodes;i++) {
        ntot += leafnodes[i].num;
        //if there are no active particles in leaf node, do nothing
        if (leafnodes[i].num == 0) continue;
        //find the near neighbours for all particles in the leaf node

#ifdef STRUCDEN
        //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
        if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) tree->FindNearestPos(leafnodes[i].cm,nnids,nnr2,opt.Nsearch);
        //otherwise distinction must be made so that only base calculation on dark matter particles
        else tree->FindNearestCheck(leafnodes[i].cm,FOFcheckpositivetype,NULL,nnids,nnr2,opt.Nsearch);
#else
        tree->FindNearestPos(leafnodes[i].cm,nnids,nnr2,opt.Nsearch);
#endif
#ifdef USEMPI
        if (opt.iLocalVelDenApproxCalcFlag==1 && NProcs > 1) {
        leafnodes[i].searchdist = sqrt(nnr2[opt.Nsearch-1]);
        //check if search region from Particle extends into other mpi domain, if so, skip particles
        bool ioverlap;
        if (opt.impiusemesh) ioverlap = (MPISearchForOverlapUsingMesh(opt,leafnodes[i].cm,leafnodes[i].searchdist)!=0);
        else ioverlap = (MPISearchForOverlap(leafnodes[i].cm,leafnodes[i].searchdist)!=0);
        if (ioverlap) continue;
        leafnodes[i].searchdist = 0;
	}
#endif
        nprocessed += leafnodes[i].num;
        for (auto j=leafnodes[i].istart;j<leafnodes[i].iend;j++)
        {
#ifdef STRUCDEN
            if (Part[j].GetType()<=0) continue;
#endif
            for (auto k=0;k<opt.Nvel;k++) {
                pqv->Push(-1, MAXVALUE);
                weight[k]=1.0;
            }
            for (auto k=0;k<opt.Nsearch;k++) {
                v2=0;
                id=nnids[k];
                if (id == j) continue;
                for (auto n=0;n<3;n++) v2+=(Part[j].GetVelocity(n)-Part[id].GetVelocity(n))*(Part[j].GetVelocity(n)-Part[id].GetVelocity(n));
                if (v2 < pqv->TopPriority()){
                    pqv->Pop();
                    pqv->Push(id, v2);
                }
            }
            Part[j].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv, weight));
        }
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] weight;
    delete pqv;
#ifdef USEOPENMP
}
#endif

#ifdef USEMPI
    //if search is fully approximative, then since particles have been localized to mpi domains in FOF groups, don't search neighbour mpi domains
    if (NProcs >1 && opt.iLocalVelDenApproxCalcFlag==1) {
    if (opt.iverbose) {
        cout<<ThisTask<<" finished local calculation in "<<MyElapsedTime(time2)<<endl;
        cout<<" fraction that is local"<<nprocessed/(float)ntot<<endl;
    }
    time2=MyGetTime();
    maxrdist=new Double_t[nbodies];
    //double maxmaxr = 0;
    for (auto i=0;i<numleafnodes;i++)
    {
        for (auto j=leafnodes[i].istart;j<leafnodes[i].iend;j++) maxrdist[j]=leafnodes[i].searchdist;
    }

    //determines export AND import numbers
    if (opt.impiusemesh) MPIGetNNExportNumUsingMesh(opt, nbodies, Part, maxrdist);
    else MPIGetNNExportNum(nbodies, Part, maxrdist);

    NNDataIn = NNDataGet = NULL;
    if (NExport>0) NNDataIn = new nndata_in[NExport];
    if (NImport>0) NNDataGet = new nndata_in[NImport];
    //build the exported particle list using NNData structures
    if (opt.impiusemesh) MPIBuildParticleNNExportListUsingMesh(opt, nbodies, Part, maxrdist);
    else MPIBuildParticleNNExportList(nbodies, Part, maxrdist);
    delete[] maxrdist;
    MPIGetNNImportNum(nbodies, tree, Part,(!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)));
    PartDataIn = PartDataGet = NULL;
    if (NExport>0) PartDataIn = new Particle[NExport];
    if (NImport>0) PartDataGet = new Particle[NImport];
    //run search on exported particles and determine which local particles need to be exported back (or imported)
    nimport=MPIBuildParticleNNImportList(opt, nbodies, tree, Part, (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)));
    int nimportsearch=opt.Nsearch;
    if (nimportsearch>nimport) nimportsearch=nimport;
    if (opt.iverbose) cout<<ThisTask<<" Searching particles in other domains "<<nimport<<endl;

    //now with imported particle list and local particle list can run proper NN search
    //first build neighbouring tree
    KDTree *treeneighbours=NULL;
    if (nimport>0) {
        treeneighbours=new KDTree(PartDataGet,nimport,1,tree->TPHYS,tree->KEPAN,100,0,0,0,period);
        treeneighbours->SetResetOrder(false);
        nprocessed=0;
    }

    //get memory useage
    GetMemUsage(opt, __func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=1));

    if (nimport>0) {
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(id,v2,nnids,nnr2,nnidsneighbours,nnr2neighbours,weight,pqx,pqv,Pval,pid2)
{
#endif

    nnids=new Int_t[opt.Nsearch];
    nnr2=new Double_t[opt.Nsearch];
    nnidsneighbours=new Int_t[opt.Nsearch];
    nnr2neighbours=new Double_t[opt.Nsearch];
    weight=new Double_t[opt.Nvel];
    pqx=new PriorityQueue(opt.Nsearch);
    pqv=new PriorityQueue(opt.Nvel);
#ifdef USEOPENMP
#pragma omp for schedule(dynamic) \
reduction(+:nprocessed)
#endif
    for (auto i=0;i<numleafnodes;i++) {
        //if there are no active particles in leaf node, do nothing
        if (leafnodes[i].num == 0) continue;
        if (leafnodes[i].searchdist == 0) continue;
        nprocessed += leafnodes[i].num;
        //find the near neighbours for all particles in the leaf node
#ifdef STRUCDEN
        //if not searching all particles in FOF then also doing baryon search then just find nearest neighbours
        if (!(opt.iBaryonSearch>=1 && opt.partsearchtype==PSTALL)) tree->FindNearestPos(leafnodes[i].cm,nnids,nnr2,opt.Nsearch);
        //otherwise distinction must be made so that only base calculation on dark matter particles
        else tree->FindNearestCheck(leafnodes[i].cm,FOFcheckpositivetype,NULL,nnids,nnr2,opt.Nsearch);
#else
        tree->FindNearestPos(leafnodes[i].cm,nnids,nnr2,opt.Nsearch);
#endif
        //fill priority queue with local particles
        for (auto j = 0; j < opt.Nsearch; j++) pqx->Push(-1, MAXVALUE);
        for (auto j = 0; j < opt.Nsearch; j++) {
            if (nnr2[j] < pqx->TopPriority()){
                pqx->Pop();
                pqx->Push(nnids[j], nnr2[j]);
            }
        }
        //search neighbouring domain and update priority queue
        treeneighbours->FindNearestPos(leafnodes[i].cm,nnidsneighbours,nnr2neighbours,nimportsearch);
        for (auto j = 0; j < nimportsearch; j++) {
            if (nnr2neighbours[j] < pqx->TopPriority()){
                pqx->Pop();
                pqx->Push(nnidsneighbours[j]+nbodies, nnr2neighbours[j]);
            }
        }
        for (auto j = 0; j < opt.Nsearch; j++) {
            nnids[j] = pqx->TopQueue();
            nnr2[j] = pqx->TopPriority();
            pqx->Pop();
        }
        for (auto j = leafnodes[i].istart; j < leafnodes[i].iend; j++)
        {
#ifdef STRUCDEN
            if (Part[j].GetType()<=0) continue;
#endif
            for (auto k=0;k<opt.Nvel;k++) {
                pqv->Push(-1, MAXVALUE);
                weight[k]=1.0;
            }
            for (auto k=0;k<opt.Nsearch;k++) {
                if (nnids[k]<nbodies) {
                    pid2=nnids[k];
                    if (pid2 == j) continue;
                    Pval = &Part[pid2];
                }
                else {
                    pid2=nnids[k]-nbodies;
                    Pval = &PartDataGet[pid2];
                }
                v2=0;for (auto n=0;n<3;n++) v2+=(Part[j].GetVelocity(n)-Pval->GetVelocity(n))*(Part[j].GetVelocity(n)-Pval->GetVelocity(n));
                if (v2 < pqv->TopPriority()){
                    pqv->Pop();
                    pqv->Push(nnids[k], v2);
                }
            }
            Part[j].SetDensity(tree->CalcSmoothLocalValue(opt.Nvel, pqv, weight));
        }
    }
    delete[] nnids;
    delete[] nnr2;
    delete[] nnidsneighbours;
    delete[] nnr2neighbours;
    delete[] weight;
    delete pqx;
    delete pqv;
#ifdef USEOPENMP
}
#endif
        delete treeneighbours;
    }
    delete[] PartDataIn;
    delete[] PartDataGet;
    delete[] NNDataIn;
    delete[] NNDataGet;
    if(opt.iverbose) {
        cout<<ThisTask<<" finished other domain search "<<MyElapsedTime(time2)<<endl;
        cout<<ThisTask<<" mpi processed fraction "<<nprocessed/(float)ntot<<endl;
    }
    }
#endif

    //free memory
    if (itreeflag) delete tree;
    if (period!=NULL) delete[] period;
}
