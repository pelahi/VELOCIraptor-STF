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
    Int_t i,j,k;
    int nthreads;
    int tid,id,pid,pid2,itreeflag=0;
    Double_t v2;
    Double_t time1,time2;
#ifndef USEMPI
    int ThisTask=0, NProcs=1;
#endif
    ///\todo alter period so arbitrary dimensions
    Double_t *period=NULL;
    if (opt.p>0) {
        period=new Double_t[3];
        for (int j=0;j<3;j++) period[j]=opt.p;
    }
    time1=MyGetTime();
    cout<<ThisTask<<": Get local velocity density"<<endl;
    //only build tree if necessary
    if (tree==NULL) {
        itreeflag=1;
        if (opt.iverbose) cout<<"Building Tree in (x) space to get local velocity density"<<endl;
        tree=new KDTree(Part,nbodies,opt.Bsize,tree->TPHYS,tree->KEPAN,1000,0,0,0,period);
    }
    if (opt.iverbose) {
        cout<<ThisTask<<" "<<"Using the following parameters to calculate velocity density using sph kernel: ";
        cout<<ThisTask<<" "<<"(Nse,Nv)="<<opt.Nsearch<<","<<opt.Nvel<<endl;
        cout<<ThisTask<<" "<<"Get velocity density using a subset of nearby physical or phase-space neighbours"<<endl;
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

    time2=MyGetTime();
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
        if (!(opt.iBaryonSearch==1 && opt.partsearchtype==PSTALL)) tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
        //otherwise distinction must be made so that only base calculation on dark matter particles
        else tree->FindNearestCriterion(i,FOFPositivetypes,NULL,nnids,nnr2,opt.Nsearch);
#else
        tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
#endif
        //once NN set is found, store maxrdist and see if particle's search radius overlaps with another mpi domain
        maxrdist[i]=sqrt(nnr2[opt.Nsearch-1]);
#ifdef SWIFTINTERFACE
        if (MPISearchForOverlapUsingMesh(libvelociraptorOpt,Part[i],maxrdist[i])==0) {
#else
        if (MPISearchForOverlap(Part[i],maxrdist[i])==0) {
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
    if (opt.iverbose) cout<<ThisTask<<" finished local calculation in "<<MyGetTime()-time2<<endl;
    time2=MyGetTime();

    //determines export AND import numbers
#ifdef SWIFTINTERFACE
    MPIGetNNExportNumUsingMesh(libvelociraptorOpt, nbodies, Part, maxrdist);
#else
    MPIGetNNExportNum(nbodies, Part, maxrdist);
#endif
    NNDataIn = new nndata_in[NExport];
    NNDataGet = new nndata_in[NImport];
    //build the exported particle list using NNData structures
#ifdef SWIFTINTERFACE
    MPIBuildParticleNNExportListUsingMesh(libvelociraptorOpt, nbodies, Part, maxrdist);
#else
    MPIBuildParticleNNExportList(nbodies, Part, maxrdist);
#endif
    MPIGetNNImportNum(nbodies, tree, Part);
    PartDataIn = new Particle[NExport];
    PartDataGet = new Particle[NImport];
    MPI_Barrier(MPI_COMM_WORLD);
    //run search on exported particles and determine which local particles need to be exported back (or imported)
    nimport=MPIBuildParticleNNImportList(nbodies, tree, Part, (!(opt.iBaryonSearch==1 && opt.partsearchtype==PSTALL)));
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
            if (!(opt.iBaryonSearch==1 && opt.partsearchtype==PSTALL)) tree->FindNearest(i,nnids,nnr2,opt.Nsearch);
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
                /*if (!(opt.iBaryonSearch==1 && opt.partsearchtype==PSTALL)) {
                    Coordinate x(Part[i].GetPosition());
                    treeneighbours->FindNearestPos(x,nnidsneighbours,nnr2neighbours,nimportsearch);
                }
                else treeneighbours->FindNearestCriterion(Part[i],FOFPositivetypes,NULL,nnidsneighbours,nnr2neighbours,nimportsearch);*/
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
    if(opt.iverbose) cout<<ThisTask<<" finished other domain search "<<MyGetTime()-time2<<endl;
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
        if (!(opt.iBaryonSearch==1 && opt.partsearchtype==PSTALL)) Part[i].SetDensity(tree->CalcVelDensityParticle(i,opt.Nvel,opt.Nsearch,1,pqx[tid],pqv[tid],&nnids[tid*opt.Nsearch],&nnr2[tid*opt.Nsearch]));
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
    cout<<ThisTask<<": finished calculation in "<<MyGetTime()-time1<<endl;
}
