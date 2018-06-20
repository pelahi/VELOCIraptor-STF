/*! \file haloproperties.cxx
 *  \brief this file contains routines to characterize the bulk properties of the halo IFF one needs to scale linking lengths

    \todo This is still not MPI compatible as it loops over local mpi thread data
 */

#include "stf.h"

/*! The FoF linking terms are adjusted by determining the physical extent of the halo and Rvir, Vcirc(Rvir), Mvir, Renc where Menclosed is [20%, 50%, 80%] of mass \n
    Units should be V=100 km/s, L=kpc, M=2.32e9 Msun so that G=1 and rhoc=3*Ho^2/8piG which gives 1.19e-7. 
    These quantites, such as average spacing and ciricular velocity at virial radius, are stored in \ref Options.ellxscale and \ref Options.ellvscale.
*/
void ScaleLinkingLengths(Options &opt,const Int_t nbodies, Particle *Part,Coordinate &cm,Coordinate &cmvel,Double_t Mtot) {
    Double_t rhoc,MaxVcirc,rlim[3], Menc[3],Renc[3], Rvir,Mvir;
    cout<<"Determining bulk halo properties to determine physical and velocity scalings"<<endl;
    AdjusttoCM(nbodies, Part, cm, cmvel, Mtot, rlim, MaxVcirc);
    Menc[0]=0.2;Menc[1]=0.5;Menc[2]=0.8;
    rhoc=1.19e-7;
    GetVirialQuantities(nbodies, Part, Mtot, rlim, rhoc, opt.virlevel, Rvir, Mvir, 3, Menc, Renc);
    //if searching only stars or gas, alter scaling to 80% enclosed mass, otherwise use virial density.
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTSTAR) {
        rlim[1]=Renc[2];MaxVcirc=sqrt(opt.G*Menc[2]*Mtot/Renc[2]);
    }
    else{
        rlim[1]=Rvir;MaxVcirc=sqrt(opt.G*Mvir/Rvir);
    }
    opt.ellxscale=sqrt((rlim[1]-rlim[0])*(rlim[1]-rlim[0]))/pow((Double_t)nbodies,1./3.);
    opt.ellvscale=MaxVcirc;
    cout<<"new scalings for length and velocity are "<<opt.ellxscale<<" "<<opt.ellvscale<<endl;
}

///\name Routines called by ScaleLinkingLengths
//@{

/*! Using initial estimate of cm get radial extent and start iterating to get best estimate of cm and move to that frame
*/
void AdjusttoCM(const Int_t nbodies, Particle *Part, Coordinate &cm, Coordinate &cmvel, const Double_t Mtot, Double_t *rlim, Double_t &MaxVcirc, Double_t tol){
    Int_t i;
    int nthreads;
    Double_t *rmin,*rmax,*rave, *vcmax,r,vc;
    Double_t change = 1e32,ri,EncMass,cmx,cmy,cmz;
    Int_t Ninside;
    Coordinate cmold=cm;
    
#ifdef USEOPENMP
#pragma omp parallel 
{
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
}
#else
    nthreads=1;
#endif
    rmin=new Double_t[nthreads];
    rmax=new Double_t[nthreads];
    rave=new Double_t[nthreads];
    vcmax=new Double_t[nthreads];

    //iterate over initial cm estimate
    r=0;
    for (int j=0;j<3;j++) r+=(Part[0].GetPosition(j)-cm[j])*(Part[0].GetPosition(j)-cm[j]);
    r=sqrt(r);
    for (int j = 0; j < nthreads; j++) {rmax[j]=r;}
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,r)
{
#pragma omp for
#endif
    for (i=1;i<nbodies;i++)
    {
        int tid;
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        r=0;
        for (int j=0;j<3;j++) r+=(Part[i].GetPosition(j)-cm[j])*(Part[i].GetPosition(j)-cm[j]);
        r=sqrt(r);
        if (r>rmax[tid])rmax[tid]=r;
    }
#ifdef USEOPENMP
}
#endif
    rlim[2]=rmax[0];
    for (int j=1;j<nthreads;j++) if (rmax[j]>rlim[2])rlim[2]=rmax[j];
    //iterate
    ri=rlim[2];
    for (int j=0;j<3;j++) cm[j]=0;
    while (change>tol)
    {
        ri*=0.9;
        // find c/m of all particles within ri
        cmx=cmy=cmz=0.;
        EncMass=0.;
        Ninside=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,r)
{
#pragma omp for reduction(+:EncMass,Ninside,cmx,cmy,cmz)

#endif
        for (i=0;i<nbodies;i++)
        {
            Double_t x,y,z;
            x = Part[i].X() - cmold[0];
            y = Part[i].Y() - cmold[1];
            z = Part[i].Z() - cmold[2];
            if (sqrt(x*x + y*y + z*z) <= ri)
            {
                cmx += Part[i].GetMass()*Part[i].X();
                cmy += Part[i].GetMass()*Part[i].Y();
                cmz += Part[i].GetMass()*Part[i].Z();
                EncMass += Part[i].GetMass();
                Ninside++;
            }
        }
#ifdef USEOPENMP
}
#endif

        cm[0]=cmx;cm[1]=cmy;cm[2]=cmz;
        for (int j=0;j<3;j++) cm[j] /= EncMass;
        // keep making radius smaller until there's
        // less than 10% of the particles inside
        if (Ninside < 0.1 * nbodies)  break;
                   
        // check for convergence
        change=fabs((cm[0]-cmold[0])/cmold[0]);
        for (int j=1;j<3;j++)if (fabs((cm[j]-cmold[j])/cmold[j]) > change) change=fabs((cm[j]-cmold[j])/cmold[j]);
        for (int j=0;j<3;j++)cmold[j]=cm[j];
    }

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i)
{
#pragma omp for
#endif
    for (i=0;i<nbodies;i++) {
        for (int k=0;k<3;k++) Part[i].SetPosition(k,Part[i].GetPosition(k)-cm[k]);
        for (int k=0;k<3;k++) Part[i].SetVelocity(k,Part[i].GetVelocity(k)-cmvel[k]);
    }
#ifdef USEOPENMP
}
#endif

    r=Part[0].Radius();
    vc=Part[0].CircularVelocity();
    for (int j = 0; j < nthreads; j++) {rmin[j]=rmax[j]=r;rave[j]=rmin[j]*Part[0].GetMass();vcmax[j]=vc;}

#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,r,vc)
{
#pragma omp for
#endif
    for (i=1;i<nbodies;i++)
    {
        int tid;
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        r=Part[i].Radius();
        vc=Part[i].CircularVelocity();
        rave[tid]+=r*Part[i].GetMass();
        if (r<rmin[tid])rmin[tid]=r;
        if (r>rmax[tid])rmax[tid]=r;
        if (vc>vcmax[tid])vcmax[tid]=vc;
    }
#ifdef USEOPENMP
}
#endif
    rlim[0]=rmin[0];
    rlim[2]=rmax[0];
    MaxVcirc=vcmax[0];
    rlim[1]=rave[0];
    for (int j=1;j<nthreads;j++) {
        if (rmin[j]<rlim[0])rlim[0]=rmin[j];
        if (rmax[j]>rlim[2])rlim[2]=rmax[j];
        if (vcmax[j]>MaxVcirc)MaxVcirc=vcmax[j];
        rlim[1]+=rave[j];
    }

    delete[] rmin;
    delete[] rmax;
    delete[] rave;
    delete[] vcmax;

    rlim[1]/=Mtot;
    rlim[0]*=0.99;rlim[2]*=1.01;
    cout<<"System (halo) has a radial extent of (min ave max)=(";for (int j=0;j<3;j++)cout<<rlim[j]<<" ";cout<<") Lunits, and a maximum circular velocity of "<<MaxVcirc<< " Vunits."<<endl;
}

/*! Get virial quantities such as mass, radius by binnning system radial about cm
    and finding radius at which enclosed mass is some specific value or the average density is some specific value (in this case virial density)
*/
void GetVirialQuantities(const Int_t nbodies, Particle *Part, const Double_t Mtot, Double_t *rlim, const Double_t rhoc, const Double_t virlevel, Double_t &Rvir, Double_t &Mvir, Int_t nenc, Double_t *Menc, Double_t *Renc){
    Int_t nbins=pow(nbodies,1./3.),i;
    int nthreads;
    Double_t deltalogr=log10(rlim[2]/rlim[0])/(Double_t)nbins,logrmin=log10(rlim[0]);
    Double_t *Mencbin,*Mbin, **mbin, *rhoavebin;
    Double_t rhovir=rhoc*virlevel;
#ifdef USEOPENMP
#pragma omp parallel 
{
    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
}
#else
    nthreads=1;
#endif
    mbin=new Double_t*[nthreads];
    for (int j=0;j<nthreads;j++) mbin[j]=new Double_t[nbins];
    Mencbin=new Double_t[nbins];
    rhoavebin=new Double_t[nbins];
    Mbin=new Double_t[nbins];
    for (int j=0;j<nthreads;j++) for (int k=0;k<nbins;k++) mbin[j][k]=0.0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i)
{
#pragma omp for
#endif
    for (i=0;i<nbodies;i++) {
        int tid;
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        Double_t logr=log10(Part[i].Radius());
        Int_t ilogr=floor((logr-logrmin)/deltalogr);
        if (ilogr<0)ilogr=0;
        if (ilogr>=nbins) ilogr=nbins-1;
        mbin[tid][ilogr]+=Part[i].GetMass();
    }
#ifdef USEOPENMP
}
#endif
    for (int j=0;j<nbins;j++) {
        Mbin[j]=mbin[0][j];
        for (int k=1;k<nthreads;k++) Mbin[j]+=mbin[k][j];
    }
    for (int k=0;k<nthreads;k++) delete[] mbin[k];
    delete[] mbin;
    Mencbin[0]=Mbin[0];
    rhoavebin[0]=Mencbin[0]/(4.0*M_PI/3.0*pow(pow(10.0,logrmin+deltalogr),3.0));
    for (int j=1;j<nbins;j++) {
        Mencbin[j]=Mencbin[j-1]+Mbin[j];
        rhoavebin[j]=Mencbin[j]/(4.0*M_PI/3.0*pow(pow(10.0,logrmin+deltalogr*(j+1)),3.0));
    }
    //find enclosed mass radius
    Int_t itemp=0,jtemp;
    for (int j=0;j<nbins-1;j++) {
        if (Mencbin[j]/Mtot<Menc[itemp]&&Mencbin[j+1]/Mtot>Menc[itemp]){
            Renc[itemp]=pow(10.0,(Menc[itemp]-Mencbin[j]/Mtot)/(Mencbin[j+1]/Mtot-Mencbin[j]/Mtot)*deltalogr+(logrmin+deltalogr*(j+1.0)));
            if (itemp==0) jtemp=j;
            itemp++;
            if (itemp==nenc) break;
        }
    }
    //find virial radius
    Rvir=rlim[2];
    Mvir=Mtot;
    for (int j=nbins-2;j>=0;j--) 
        if (rhoavebin[j]/rhovir>1.0&&rhoavebin[j+1]/rhovir<1.0) {
            Rvir=pow(10.0,(1.0-rhoavebin[j]/rhovir)/(rhoavebin[j+1]/rhovir-rhoavebin[j]/rhovir)*deltalogr+(logrmin+deltalogr*(j+1.0)));
            Mvir=Mencbin[j]+(Mencbin[j+1]-Mencbin[j])/deltalogr*(log10(Rvir)-(logrmin+deltalogr*(j+1.0)));
            break;
        }
    delete[] Mbin;
    delete[] Mencbin;
    delete[] rhoavebin;
    cout<<"Halo has mass frac of "<<Mvir/Mtot<<" at virial radius "<<Rvir<<endl;
    cout<<"And (Mfrac,Rmfrac) of ";for (int j=0;j<nenc;j++) cout<<"("<<Menc[j]<<","<<Renc[j]<<") ";cout<<endl;
}
//@}
