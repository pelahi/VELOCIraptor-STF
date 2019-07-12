/*! \file localbgcomp.cxx
 *  \brief this file contains routines to compare the local velocity density function of a particle relative to that predicted by the backgound and then calculates the normalized deviation or outlier value ellprob
    \todo once ratio is calculated, must figure out best way to do global mpi fit. Probably best to combine the binned distribution and fit that
 */

#include "stf.h"

/*! This calculates the logarithmic ratio of the measured velocity density and the expected velocity density assuming a bg muiltivariate gaussian distribution
    \todo must adjust interpolation scheme so that if NN has cells in a neighbouring MPI domain, the information is stored locally. This may require a rewrite
    of the grid cell structure or the near neighbour list so that if grid cell has NN in another processor, actually physically store the information cm, cmvel, veldisp
    locally to that grid cell. Another option is to determine all cells that are NN of a cell in another mpi's domain, build a grid export list that contains the relevant information
    and an index that is easily accessible when searching the list of NN cells locally.
*/
void GetDenVRatio(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp)
{
    Int_t i;
    int nthreads,tid;
#ifndef USEMPI
    int ThisTask=0;
#endif
    Double_t **dist;
    Double_t norm=pow(2.0*M_PI,-1.5);
    Int_t **nn;

    Double_t w,wsum,maxdist,sv,vsv,fbg,tempdenv;
    Coordinate vp,vmweighted;
    Matrix isvweighted;
    Particle *ptemp;
    KDTree *tree;

    if (opt.iverbose>=2) cout<<ThisTask<<" Now calculate denvratios using grid"<<endl;
    //take inverse for interpolation
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<ngrid;i++) gveldisp[i]=gveldisp[i].Inverse();
#ifdef USEOPENMP
}
#endif

    //build grid tree so that one can find nearest cells for each particle
    //if using MPI since number of cells is far fewer than number of particles, simple gather collect all the data so that each processor has access to it
#ifdef USEMPI
    if(opt.iSingleHalo) {
    Ngridlocal=ngrid;
    MPI_Allreduce(&ngrid,&Ngridtotal,1,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
    mpi_grid=new GridCell[Ngridtotal];
    mpi_gvel=new Coordinate[Ngridtotal];
    mpi_gveldisp=new Matrix[Ngridtotal];
    MPIBuildGridData(ngrid, grid, gvel, gveldisp);
    delete[] grid;
    delete[] gvel;
    delete[] gveldisp;
    ngrid=Ngridtotal;
    grid=mpi_grid;
    gvel=mpi_gvel;
    gveldisp=mpi_gveldisp;
    }
#endif
    ptemp=new Particle[ngrid];
    for (i=0;i<ngrid;i++) ptemp[i]=Particle(1.0,grid[i].xm[0],grid[i].xm[1],grid[i].xm[2],0.0,0.0,0.0,i);
    tree=new KDTree(ptemp,ngrid,1,tree->TPHYS);

#ifndef USEOPENMP
    nthreads=1;
#else
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    dist=new Double_t*[nthreads];
    for (int j=0;j<nthreads;j++)dist[j]=new Double_t[MAXNGRID+1];
    nn=new Int_t*[nthreads];
    for (int j=0;j<nthreads;j++)nn[j]=new Int_t[MAXNGRID+1];

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,w,wsum,sv,vsv,fbg,vp,maxdist,vmweighted,isvweighted,tid,tempdenv)
{
#pragma omp for schedule(dynamic) nowait
#endif
    for (i=0;i<nbodies;i++)
    {
        tempdenv=Part[i].GetDensity()/opt.Nsearch;
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif

        //try inverse distance weighting scheme based using Shepard's method.
        fbg=0.;
        wsum=0.;
        maxdist=0.;
        vmweighted[0]=vmweighted[1]=vmweighted[2]=0.;
        for (int j=0;j<3;j++) for (int k=0;k<3;k++) isvweighted(j,k)=0.0;
        Coordinate xpos(Part[i].GetPosition());
        tree->FindNearestPos(xpos,nn[tid],dist[tid],MAXNGRID+1);
        for (int j=0;j<=MAXNGRID;j++) {
           dist[tid][j]=sqrt(dist[tid][j]+1e-16);
           if (dist[tid][j]>maxdist)maxdist=dist[tid][j];
        }
        for (int j=0;j<=MAXNGRID;j++) {
            w=(maxdist-dist[tid][j])/(maxdist*dist[tid][j]);w=w*w;
            //w=1.0/dist[tid][j];
            wsum+=w;
            vmweighted=vmweighted+gvel[ptemp[nn[tid][j]].GetID()]*w;
            isvweighted=isvweighted+gveldisp[ptemp[nn[tid][j]].GetID()]*w;
        }
        vmweighted=vmweighted*(1.0/wsum);
        isvweighted=isvweighted*(1.0/wsum);
        sv=sqrt(abs(isvweighted.Det()));
        for (int m=0;m<3;m++) vp[m]=Part[i].GetVelocity(m)-vmweighted[m];
        vsv=0.;for (int m=0;m<3;m++) for (int n=0;n<3;n++) vsv+=vp[m]*vp[n]*isvweighted(m,n);
        fbg=log(sv)-0.5*vsv;
        Part[i].SetPotential(log(tempdenv)-log(norm)-fbg);
    }
#ifdef USEOPENMP
}
#endif
    if (opt.iverbose>=2) cout<<ThisTask<<" Done"<<endl;
    delete[] gvel;
    delete[] gveldisp;
    delete tree;
    delete[] ptemp;
    delete[] grid;
}


void DetermineDenVRatioDistribution(Options &opt,const Int_t nbodies, Particle *Part, Double_t &meanr,Double_t &sdlow,Double_t &sdhigh, int sublevel)
{
    Int_t i,nbins,iprob,jprob;
    Double_t mtot,mtotpeak,*rbin,deltar,maxprob, minprob,rmin,rmax;
    Double_t *xbin;
    Double_t **omp_rbin;//, **omp_xbin;
    int nthreads=1,tid;
    Double_t w;
    Int_t ir;
#ifdef USEOPENMP
#pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
#endif
    //to determine initial number of bins using modified Sturges' formula
    nbins = ceil(log10((Double_t)nbodies)/log10(2.0)+1)*4;
    //rbin=new Double_t[nbins];
    //xbin=new Double_t[nbins];
    omp_rbin=new Double_t*[nthreads];
    for (i=0;i<nthreads;i++) {
        omp_rbin[i]=new Double_t[nbins];
    }

    //deterrmine average, rmin,rmax and variance about mean
    rmin=rmax=Part[0].GetPotential();
#ifdef USEOPENMP
    Double_t rmina[nthreads],rmaxa[nthreads];
    for (int j=0;j<nthreads;j++) rmina[j]=rmaxa[j]=Part[0].GetPotential();
#endif

#ifdef USEOPENMP
    if (nbodies>ompperiodnum) {
#pragma omp parallel default(shared) \
private(i,tid)
{
#pragma omp for
    for (i=1;i<nbodies;i++) {
        tid=omp_get_thread_num();
        if (rmina[tid]>Part[i].GetPotential())rmina[tid]=Part[i].GetPotential();
        if (rmaxa[tid]<Part[i].GetPotential())rmaxa[tid]=Part[i].GetPotential();
    }
}
    rmin=rmina[0];rmax=rmaxa[0];
    for (int j=1;j<nthreads;j++) {
        if (rmin>rmina[j])rmin=rmina[j];
        if (rmax<rmaxa[j])rmax=rmaxa[j];
    }
    }
    else {
    for (i=1;i<nbodies;i++) {
        if (rmin>Part[i].GetPotential())rmin=Part[i].GetPotential();
        if (rmax<Part[i].GetPotential())rmax=Part[i].GetPotential();
    }
    }
#else
    for (i=1;i<nbodies;i++) {
        if (rmin>Part[i].GetPotential())rmin=Part[i].GetPotential();
        if (rmax<Part[i].GetPotential())rmax=Part[i].GetPotential();
    }
#endif

    //now bin data and find initial estimates for most probable value and the FWHM on either side of the most probable value
    //deltar=(rmax-rmin)/(Double_t)nbins;
    deltar=(4.0*fabs(rmin))/(Double_t)nbins;
    rmin-=deltar*0.025;
    deltar*=1.05;
    //for (i=0;i<nbins;i++) rbin[i]=0;
    for (int j=0;j<nthreads;j++) for (i=0;i<nbins;i++) omp_rbin[j][i]=0;
    mtot=0;
#ifdef USEOPENMP
    if (nbodies>ompperiodnum) {
#pragma omp parallel default(shared)
{
#pragma omp for private(i,tid,w,ir) reduction(+:mtot)
    for (i=0;i<nbodies;i++) {
        tid=omp_get_thread_num();
        ir=(Int_t)((Part[i].GetPotential()-rmin)/deltar);
        //mass weighted
#ifdef NOMASSWEIGHT
        w=1.0;
#else
        w=Part[i].GetMass();
#endif
        if (ir<nbins) {
            omp_rbin[tid][ir]+=w;
            //rbin[ir]+=w;
            mtot+=w;
        }
    }
}
    }
    else {
    for (i=0;i<nbodies;i++) {
        tid=0;
        ir=(Int_t)((Part[i].GetPotential()-rmin)/deltar);
        //mass weighted
#ifdef NOMASSWEIGHT
        w=1.0;
#else
        w=Part[i].GetMass();
#endif
        if (ir<nbins) {
            omp_rbin[tid][ir]+=w;
            mtot+=w;
        }
    }
    }
#else
    for (i=0;i<nbodies;i++) {
        tid=0;
        ir=(Int_t)((Part[i].GetPotential()-rmin)/deltar);
        //mass weighted
#ifdef NOMASSWEIGHT
        w=1.0;
#else
        w=Part[i].GetMass();
#endif
        if (ir<nbins) {
            omp_rbin[tid][ir]+=w;
            mtot+=w;
        }
    }
#endif
    for (int j=1;j<nthreads;j++) for (i=0;i<nbins;i++) omp_rbin[0][i]+=omp_rbin[j][i];
    rbin=&omp_rbin[0][0];

    maxprob=0.;
    for (i=0;i<nbins;i++) {
        if (rbin[i]>maxprob) maxprob=rbin[iprob=i];
    }

    meanr=(iprob+0.5)*deltar+rmin;
    //find first estimate of sdlow by going from rmin to prob and when have some expect fraction of the probability
    Double_t sl=1.0;
    for (i=iprob;i>=0;i--) {
        if (rbin[i]<=exp(-0.5*sl*sl)*rbin[iprob]) {
            jprob=i;
            sdlow=(meanr-(((exp(-0.5*sl*sl)*rbin[iprob]-rbin[jprob])/(rbin[jprob+1]-rbin[jprob])+jprob+0.5)*deltar+rmin))/sl;
            break;
        }
        if (i==0) {
            sdlow=(iprob-jprob)*deltar/sl;
        }
    }
    for (i=iprob;i<nbins;i++){
        if (rbin[i]<=exp(-0.5*sl*sl)*rbin[iprob]) {
            jprob=i;
            sdhigh=((((exp(-0.5*sl*sl)*rbin[iprob]-rbin[jprob-1])/(rbin[jprob]-rbin[jprob-1])+jprob+0.5)*deltar+rmin)-meanr)/sl;
            break;
        }
        if (i==nbins-1) {
            jprob=nbins-1;
            sdhigh=(jprob-iprob)*deltar/sl;
        }
    }
    //if object is small or bg search (ie sublevel==-1, then to keep statistics high, use preliminary determination of the variance and mean.
    if (nbodies<2*MINSUBSIZE) {
        if (opt.iverbose>=2) printf("Using meanr=%e sdlow=%e sdhigh=%e\n",meanr,sdlow,sdhigh);
        return;
    }
    //now rebin around most probable over sl in either direction to be used to estimate dispersion
    //and gradually increase region till region encompases over 50% of the mass or particle numbers
    GMatrix W(nbins,nbins);
    rbin=new Double_t[nbins];
    do {
        mtotpeak=0;
        rmin=(meanr-sl*sdlow);
        rmax=(meanr+sl*sdhigh);
        int npeak=0;
        for (i=0;i<nbodies;i++) if (Part[i].GetPotential()>=rmin&&Part[i].GetPotential()<rmax) npeak++;
        //delete[] rbin;
        //delete[] xbin;
        //for (i=0;i<nthreads;i++) delete[] omp_rbin[i];
        //once have initial estimates of variance bin using Scott's formula
        //deltar=3.5*sdlow/pow(nbodies,1./3.);
        deltar=3.5*sqrt(sdlow*sdlow+sdhigh*sdhigh)/pow(npeak,1./3.);
        nbins=ceil((rmax-rmin)/deltar+1);
        W=GMatrix(nbins,nbins);
		delete[] rbin;
        rbin=new Double_t[nbins];
        //for (i=0;i<nthreads;i++) omp_rbin[i]=new Double_t[nbins];
        //xbin=new Double_t[nbins];
        for (i=0;i<nbins;i++)rbin[i]=0;
        //for (int j=0;j<nthreads;j++) for (i=0;i<nbins;i++) omp_rbin[j][i]=0;
        for (int j=0;j<nbins;j++) for (int k=0;k<nbins;k++) W(j,k)=0.;
/*#ifdef USEOPENMP
#pragma omp parallel default(shared)
{
#pragma omp for private(i) reduction(+ : mtotpeak)
#endif
*/
        for (i=0;i<nbodies;i++) {
//            int tid=0;
//#ifdef USEOPENMP
//            tid=omp_get_thread_num();
//#endif
            if (Part[i].GetPotential()>=rmin&&Part[i].GetPotential()<rmax) {
                ir=(Int_t)((Part[i].GetPotential()-rmin)/deltar);
#ifdef NOMASSWEIGHT
                w=1.0;
#else
                w=Part[i].GetMass();
#endif
                rbin[ir]+=w;
                W(ir,ir)+=w*w;
                //omp_rbin[tid][ir]+=w;
                mtotpeak+=w;
            }
        }
/*#ifdef USEOPENMP
}
#endif*/
        sl*=1.25;
    }while (mtotpeak/mtot<0.2);
    GMatrix covar(nbins,nbins);
    //add bins together
    //for (int j=1;j<nthreads;j++) for (i=0;i<nbins;i++) omp_rbin[0][i]+=omp_rbin[j][i];
    //rbin=omp_rbin[0];
    xbin=new Double_t[nbins];
    maxprob=0.;
    minprob=MAXVALUE;
    for (i=0;i<nbins;i++) {
        if (rbin[i]>maxprob)maxprob=rbin[iprob=i];
        //if (rbin[i]<minprob&&rbin[i]>0.)minprob=rbin[i];
        if (W(i,i)<minprob&&rbin[i]>0.) minprob=W(i,i);
        xbin[i]=(i+0.5)*deltar+rmin;
    }
    for (i=0;i<nbins;i++) {
        if (rbin[i]!=0) W(i,i)=1.0/W(i,i);
        else W(i,i)=1.0/minprob;//simple Poisson weighting with value in bin as mean and Poisson distributed about this mean
    }
    meanr=(iprob+0.5)*deltar+rmin;
    sl=0.9;
    for (i=iprob;i>=0;i--) {
        if (rbin[i]<=exp(-0.5*sl*sl)*rbin[iprob]) {
            jprob=i;
            sdlow=(meanr-(((exp(-0.5*sl*sl)*rbin[iprob]-rbin[jprob])/(rbin[jprob+1]-rbin[jprob])+jprob+0.5)*deltar+rmin))/sl;
            break;
        }
        /*
        if (i==0) {
            sdlow=iprob*deltar/sl;
        }
        */
    }
    for (i=iprob;i<nbins;i++){
        if (rbin[i]<=exp(-0.5*sl*sl)*rbin[iprob]) {
            jprob=i;
            sdhigh=((((exp(-0.5*sl*sl)*rbin[iprob]-rbin[jprob-1])/(rbin[jprob]-rbin[jprob-1])+jprob+0.5)*deltar+rmin)-meanr)/sl;
            break;
        }
        if (i==nbins-1) {
            jprob=nbins-1;
            sdhigh=(jprob-iprob)*deltar/sl;
        }
    }
    //adjust sdhigh to sdlow due to assymetry
    sdhigh=sdlow;
    //again, if number of particles is low (and so bin statisitics is poor) use initial estimate
    if (nbodies<16*MINSUBSIZE||sublevel==-1) {
        if (opt.iverbose>=2) printf("Using meanr=%e sdlow=%e sdhigh=%e\n",meanr,sdlow,sdhigh);
        return;
    }

    //check if mean consistent with zero, then set to zero for fitting
    //if ((meanr-sdlow)*(meanr+sdhigh)<0) meanr=0;

    //now have initial estimates of paramters, try nonlinear ls fit to data below prob and above
    Int_t iflag,iit=0,nparams=4;
    Double_t chi2,oldchi2;
    Double_t *params=new Double_t[nparams];
    //five sets of fix parameter choices so that get optimal fit given bad data.
    int **fixp,nfix,itemp;
    int nfits=8;
    fixp= new int*[nfits];
    for (int i=0;i<8;i++) fixp[i]=new int[nparams];
    struct math_function fitfunc,*difffuncs;
    difffuncs=new math_function[nparams];

    nparams=4;
    fitfunc.function=SkewGauss;
    difffuncs[0].function=DiffSkewGaussAmp;
    difffuncs[1].function=DiffSkewGaussMean;
    difffuncs[2].function=DiffSkewGaussVar;
    difffuncs[3].function=DiffSkewGaussSkew;
    params[0]=maxprob;
    params[1]=meanr;
    params[2]=sdhigh*sdhigh*0.8;//assume conservative dispersion
    params[3]=1.0;

    itemp=0;
    fixp[itemp][0]=1;fixp[itemp][1]=1;fixp[itemp][2]=0;fixp[itemp][3]=1;itemp++;
    fixp[itemp][0]=1;fixp[itemp][1]=0;fixp[itemp][2]=0;fixp[itemp][3]=1;itemp++;
    fixp[itemp][0]=0;fixp[itemp][1]=0;fixp[itemp][2]=0;fixp[itemp][3]=1;itemp++;
    fixp[itemp][0]=1;fixp[itemp][1]=1;fixp[itemp][2]=1;fixp[itemp][3]=0;itemp++;
    fixp[itemp][0]=0;fixp[itemp][1]=1;fixp[itemp][2]=1;fixp[itemp][3]=1;itemp++;
    fixp[itemp][0]=1;fixp[itemp][1]=0;fixp[itemp][2]=0;fixp[itemp][3]=1;itemp++;
    fixp[itemp][0]=1;fixp[itemp][1]=0;fixp[itemp][2]=0;fixp[itemp][3]=0;itemp++;
    fixp[itemp][0]=0;fixp[itemp][1]=0;fixp[itemp][2]=0;fixp[itemp][3]=0;itemp++;

    if (opt.iverbose>=2) printf("Initial estimate: mu=%e var=%e \n",params[1],sqrt(params[2]));
    nfits=8;
    oldchi2=MAXVALUE;
    for (int i=0;i<nfits;i++) {
        chi2=FitNonLinLS(fitfunc, difffuncs, nparams, params, covar, nbins, xbin, rbin, &W,  1e-2, 0.95, fixp[i],1,20);
        int ifitfail=0;
        for (int j=0;j<nparams;j++) ifitfail+=std::isnan(params[j]);
        ifitfail+=(params[2]<=0);
        ifitfail+=(params[3]<=0);
        if (chi2<oldchi2&&chi2!=-1&&std::isnan(chi2)==0&&ifitfail==0) {
            meanr=params[1];sdlow=sqrt(params[2]*params[3]);sdhigh=sqrt(params[2]);
            nfix=0;for (int j=0;j<nparams;j++) nfix+=(fixp[i][j]==1);
            oldchi2=chi2;
            if(opt.iverbose>2) printf("chi2/dof=%e/%d, A=%e mu=%e var=%e s=%e\n",chi2,nbins-(nparams-nfix)-1,params[0],params[1],sqrt(params[2]),sqrt(params[3]));
        }
        else if (oldchi2<chi2) break;
        else {
            if (opt.iverbose>2)printf("fit failed, using previous values\n");
            params[0]=maxprob;params[1]=meanr;params[2]=sdhigh*sdhigh;params[3]=(sdlow*sdlow)/(sdhigh*sdhigh);
        }
    }

    if (opt.iverbose>=2) printf("Using meanr=%e sdlow=%e sdhigh=%e\n",meanr,sdlow,sdhigh);
    //free memory
    delete[] xbin;
    delete[] rbin;
    //rbin=NULL;
    for (int j=0;j<nthreads;j++) delete[] omp_rbin[j];
    delete[] omp_rbin;
}

/*! Calculates the normalized deviations from the mean of the dominated population.
    \todo note that before had FOFSTPROB set density to probability, but here set to ell, the normalized logaritmic "distance" from predicted maxwellian velocity density)
    but could add routine that transforms these values to probablity if necessary.

*/
Int_t GetOutliersValues(Options &opt, const Int_t nbodies, Particle *Part, int sublevel)
{
    Int_t i,nsubset;
    int nthreads;
    Double_t temp, mtot=0.0;
#ifndef USEMPI
    int ThisTask=0;
#endif
    if (opt.iverbose>=2) cout<<ThisTask<<" Now get average in grid cell and find outliers"<<endl;
    //printf("Using GLOBAL values to characterize the distribution and determine the normalized values used to determine outlier likelihood\n");
    Double_t globalave,globalvar,globalmostprob,globalsdlow,globalsdhigh;

    DetermineDenVRatioDistribution(opt,nbodies,Part,globalmostprob,globalsdlow,globalsdhigh, sublevel);

    Double_t temp1,temp2,temp3, tempell;
    temp1=globalmostprob;
    temp2=1.0/(globalsdhigh);
    temp3=1.0/(globalsdlow);
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,tempell) \
firstprivate(temp1,temp2,temp3)
{
#pragma omp for reduction(+:nsubset)
#endif
    for (i=0;i<nbodies;i++) {
        tempell=(Part[i].GetPotential()-globalmostprob);
        if (tempell>0) tempell*=temp2;
        else tempell*=temp3;
        Part[i].SetPotential(tempell);
        nsubset+=(Part[i].GetPotential()>opt.ellthreshold);
    }
#ifdef USEOPENMP
}
#endif
    if (opt.iverbose>=2) cout<<ThisTask<<" Done"<<endl;
    return nsubset;
}
