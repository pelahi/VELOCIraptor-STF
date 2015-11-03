/*! \file properties.cxx
 *  \brief this file contains routines to calculate bulk properties and sort list according to them
 */

#include "baryoniccontent.h"

double MyGetTime(){
#ifdef USEOPENMP
    return omp_get_wtime();
#else
    return (clock() /(60.* (double)CLOCKS_PER_SEC));
#endif
}

int Double_t_cmp(const void *a, const void *b)
{
        Double_t aa = (*(Double_t*)a);
        Double_t bb = (*(Double_t*)b);
        if (aa > bb) return 1;
        else if (aa < bb) return -1;
        else return 0;
}
///get the stats of data
void GetStats(Int_t ndata,Double_t *data,Double_t *results,int nquants, Double_t *quants){
    Int_t i;
    Double_t mean,std,mquant,iquant,gquant;
    mean=std=0;
    qsort(data,ndata,sizeof(Double_t),Double_t_cmp);
    if (ndata<=omppropnum){
        for (i=0;i<ndata;i++) {
            mean+=data[i];
            std+=data[i]*data[i];
        }
        mean/=(Double_t)ndata;
        std/=(Double_t)ndata;
        std-=mean*mean;
        results[0]=mean;
        results[1]=sqrt(std);
        for (i=0;i<nquants;i++){
            mquant=0.4+quants[i]*(1-0.4-0.4);
            iquant = floor(ndata*quants[i]+mquant);
            gquant=ndata*quants[i]+mquant-iquant;
            if (iquant<=ndata-1 &&iquant>=0)
                results[i+2]=(1.0-gquant)*data[(int)iquant]+gquant*data[(int)iquant];
            else if (iquant<0)
                results[i+2]=data[0];
            else if (iquant>ndata-1)
                results[i+2]=data[ndata-1];
        }
    }
    else {
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i)
{
    #pragma omp for reduction(+:mean,std)
#endif
        for (i=0;i<ndata;i++) {
            mean+=data[i];
            std+=data[i]*data[i];
        }
#ifdef USEOPENMP
}
#endif
        mean/=(Double_t)ndata;
        std/=(Double_t)ndata;
        std-=mean*mean;
        results[0]=mean;
        results[1]=sqrt(std);
        for (i=0;i<nquants;i++){
            mquant=0.4+quants[i]*(1-0.4-0.4);
            iquant = floor(ndata*quants[i]+mquant);
            gquant=ndata*quants[i]+mquant-iquant;
            if (iquant<=ndata-1&&iquant>=0)
                results[i+2]=(1.0-gquant)*data[(int)iquant]+gquant*data[(int)iquant];
            else if (iquant<0) results[i+2]=data[0];
            else if (iquant>ndata-1) results[i+2]=data[ndata-1];
        }
    }
}

void AdjustForPeriod(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp){
    Int_t i,j,k;
    for (i=0;i<ngroup;i++) {
        for (j=1;j<hp[i].NumberofParticles;j++) {
            for (k=0;k<3;k++) {
                if (Part[j+hp[i].noffset].GetPosition(k)-Part[hp[i].noffset].GetPosition(k)>0.5*opt.p)
                    Part[j+hp[i].noffset].SetPosition(k,Part[j+hp[i].noffset].GetPosition(k)-opt.p);
                else if (Part[j+hp[i].noffset].GetPosition(k)-Part[hp[i].noffset].GetPosition(k)<-0.5*opt.p)
                    Part[j+hp[i].noffset].SetPosition(k,Part[j+hp[i].noffset].GetPosition(k)+opt.p);
            }
        }
    }
}

///Gets cm
///can be used to cm of a particular type of particle, where assumes that particles in a structure are sorted according to type and then binding energy. Just need to alter noffset and numingroup appropriately
void GetCM(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp, int ptype)
{
    cout<<"Get CM"<<endl;
    Particle *Pval;
    Int_t i,j,k,noffset;
    Double_t ri,rcmv,r2,cmx,cmy,cmz,EncMass,Ninside;
    Double_t vc,rc,x,y,z,vx,vy,vz;
    Coordinate cmold(0.),cmref;
    Double_t change=MAXVALUE,tol=1e-2;
    Int_t ii,icmv;
    Double_t virval=log(opt.virlevel*opt.rhobg);
    Int_t *numingroup=new Int_t[ngroup];
    Double_t t1=MyGetTime();

    Coordinate *blahcms=new Coordinate[ngroup];
    Coordinate *blahcmvels=new Coordinate[ngroup];
    for (i=0;i<ngroup;i++) blahcms[i]=pdata[i].gcm;
    for (i=0;i<ngroup;i++) blahcmvels[i]=pdata[i].gcmvel;

#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) {
        if (ptype==-1) numingroup[i]=hp[i].NumberofParticles;
        else numingroup[i]=hp[i].NumofType[ptype];
        pdata[i].num=numingroup[i];
    }
#ifdef USEOPENMP
}
#endif
    //for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,Pval,ri,rcmv,r2,cmx,cmy,cmz,EncMass,Ninside,cmold,change,tol,x,y,z,vc,rc,vx,vy,vz,noffset)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (numingroup[i]<=omppropnum&&numingroup[i]>0)
    {
        //set offset correctly
        noffset=hp[i].noffset;
        for (j=1;j<=ptype;j++) noffset+=hp[i].AllNumofType[j-1];
        for (k=0;k<3;k++) pdata[i].gcm[k]=pdata[i].gcmvel[k]=0;
        pdata[i].gmass=0.0;
        EncMass=0.;
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            EncMass+=(*Pval).GetMass();
            for (k=0;k<3;k++) {
                pdata[i].gcm[k]+=(*Pval).GetPosition(k)*(*Pval).GetMass();
                pdata[i].gcmvel[k]+=(*Pval).GetVelocity(k)*(*Pval).GetMass();
            }
        }
        pdata[i].gmass=EncMass;
        for (k=0;k<3;k++){pdata[i].gcm[k]*=(1.0/pdata[i].gmass);pdata[i].gcmvel[k]*=(1.0/pdata[i].gmass);}
        pdata[i].gsize=0;
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            r2=0.0;
            for (k=0;k<3;k++) r2+=(pdata[i].gcm[k]-(*Pval).GetPosition(k))*(pdata[i].gcm[k]-(*Pval).GetPosition(k));
            if (sqrt(r2)>pdata[i].gsize)pdata[i].gsize=sqrt(r2);
        }
        //iterate for better cm if group large enough
        //iteration is done by slowly shrinking the radial sphere used to estimate the cm
        cmold[0]=cmold[1]=cmold[2]=0.;
        change=MAXVALUE;tol=1e-2;
        if (numingroup[i]>=50) {
            ri=pdata[i].gsize;
            ri=ri*ri;
            cmold=pdata[i].gcm;
            rcmv=ri;
            while (true)
            {
                ri*=opt.cmadjustfac;
                // find c/m of all particles within ri
                cmx=cmy=cmz=0.;
                EncMass=0.;
                Ninside=0;
                for (j=0;j<numingroup[i];j++)
                {
                    Pval=&Part[j+noffset];
                    x = (*Pval).X() - cmold[0];
                    y = (*Pval).Y() - cmold[1];
                    z = (*Pval).Z() - cmold[2];
                    if ((x*x + y*y + z*z) <= ri)
                    {
                        cmx += (*Pval).GetMass()*(*Pval).X();
                        cmy += (*Pval).GetMass()*(*Pval).Y();
                        cmz += (*Pval).GetMass()*(*Pval).Z();
                        EncMass += (*Pval).GetMass();
                        Ninside++;
                    }
                }
                if (Ninside > opt.cmfrac * numingroup[i] && EncMass>0.) {
                    pdata[i].gcm[0]=cmx;pdata[i].gcm[1]=cmy;pdata[i].gcm[2]=cmz;
                    for (k=0;k<3;k++) pdata[i].gcm[k] /= EncMass;
                    cmold=pdata[i].gcm;
                    rcmv=ri;
                }
		        else break;
            }
            cmx=cmy=cmz=EncMass=0.;
            for (j=0;j<numingroup[i];j++)
            {
                Pval=&Part[j+noffset];
                x = (*Pval).X() - pdata[i].gcm[0];
                y = (*Pval).Y() - pdata[i].gcm[1];
                z = (*Pval).Z() - pdata[i].gcm[2];
                if ((x*x + y*y + z*z) <= rcmv)
                {
                    cmx += (*Pval).GetMass()*(*Pval).Vx();
                    cmy += (*Pval).GetMass()*(*Pval).Vy();
                    cmz += (*Pval).GetMass()*(*Pval).Vz();
                    EncMass += (*Pval).GetMass();
                }
            }
            pdata[i].gcmvel[0]=cmx;pdata[i].gcmvel[1]=cmy;pdata[i].gcmvel[2]=cmz;
            for (k=0;k<3;k++) pdata[i].gcmvel[k] /= EncMass;
        }
    }
#ifdef USEOPENMP
}
#endif

    for (i=0;i<ngroup;i++) if (numingroup[i]>omppropnum)
    {
        //set offset correctly
        noffset=hp[i].noffset;
        for (j=1;j<=ptype;j++) noffset+=hp[i].AllNumofType[j-1];

        for (k=0;k<3;k++) pdata[i].gcm[k]=pdata[i].gcmvel[k]=0;
        pdata[i].gmass=pdata[i].gmaxvel=0.0;
        EncMass=cmx=cmy=cmz=0.;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,Pval)
{
    #pragma omp for reduction(+:EncMass,cmx,cmy,cmz)
#endif
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            EncMass+=(*Pval).GetMass();
            cmx+=(*Pval).X()*(*Pval).GetMass();
            cmy+=(*Pval).Y()*(*Pval).GetMass();
            cmz+=(*Pval).Z()*(*Pval).GetMass();
        }
#ifdef USEOPENMP
}
#endif
        pdata[i].gcm[0]=cmx;pdata[i].gcm[1]=cmy;pdata[i].gcm[2]=cmz;
        pdata[i].gmass=EncMass;
        for (k=0;k<3;k++){pdata[i].gcm[k]*=(1.0/pdata[i].gmass);pdata[i].gcmvel[k]*=(1.0/pdata[i].gmass);}
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            for (k=0;k<3;k++) {
                Pval->SetPosition(k,(*Pval).GetPosition(k)-pdata[i].gcm[k]);
            }
        }
        qsort(&Part[noffset], numingroup[i], sizeof(Particle), RadCompare);
        ri=Part[noffset+numingroup[i]-1].Radius();
        ri=ri*ri;
        //iterate for better cm if group large enough
        cmold[0]=cmold[1]=cmold[2]=0.;
        change=MAXVALUE;tol=1e-2;
        cmref=pdata[i].gcm;//cmold=pdata[i].gcm;
        rcmv=ri;
        ii=numingroup[i];
        while (true)
        {
            ri*=opt.cmadjustfac;
            //ii*=opt.cmadjustfac;
            // find c/m of all particles within ri
            cmx=cmy=cmz=0.;
            EncMass=0.;
            Ninside=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(j,Pval,x,y,z)
{
#pragma omp for reduction(+:EncMass,Ninside,cmx,cmy,cmz)
#endif
            for (j=0;j<numingroup[i];j++)
            {
                Pval=&Part[j+noffset];
                x = (*Pval).X() - cmold[0];
                y = (*Pval).Y() - cmold[1];
                z = (*Pval).Z() - cmold[2];
                if ((x*x + y*y + z*z) <= ri)
                {
                    cmx += (*Pval).GetMass()*(*Pval).X();
                    cmy += (*Pval).GetMass()*(*Pval).Y();
                    cmz += (*Pval).GetMass()*(*Pval).Z();
                    EncMass += (*Pval).GetMass();
                    Ninside++;
                }
            }
#ifdef USEOPENMP
}
#endif
            x = Part[noffset+ii-1].X() - cmold[0];
            y = Part[noffset+ii-1].Y() - cmold[1];
            z = Part[noffset+ii-1].Z() - cmold[2];
            //ri=x*x+y*y+z*z;
            if (Ninside > opt.cmfrac * numingroup[i] && EncMass>0.) {
                cmold[0]=cmx;cmold[1]=cmy;cmold[2]=cmz;
                for (k=0;k<3;k++) cmold[k] /= EncMass;
                rcmv=ri;
                //icmv=ii;
            }
            else break;
        }
        double oldenc=EncMass;
        double oldninside=Ninside;
        double oldri=ri;
        cmx=cmy=cmz=EncMass=0.;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(j,Pval,x,y,z)
{
    #pragma omp for reduction(+:EncMass,cmx,cmy,cmz)
#endif
        for (j=0;j<numingroup[i];j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - cmold[0];
            y = (*Pval).Y() - cmold[1];
            z = (*Pval).Z() - cmold[2];
            if ((x*x + y*y + z*z) <= rcmv)
            {
                cmx += (*Pval).GetMass()*(*Pval).Vx();
                cmy += (*Pval).GetMass()*(*Pval).Vy();
                cmz += (*Pval).GetMass()*(*Pval).Vz();
                EncMass += (*Pval).GetMass();
            }
            x = (*Pval).X() + pdata[i].gcm[0];
            y = (*Pval).Y() + pdata[i].gcm[1];
            z = (*Pval).Z() + pdata[i].gcm[2];
            Pval->SetPosition(x,y,z);
        }
#ifdef USEOPENMP
}
#endif
        for (k=0;k<3;k++) pdata[i].gcm[k]+=cmold[k];
        pdata[i].gcmvel[0]=cmx;pdata[i].gcmvel[1]=cmy;pdata[i].gcmvel[2]=cmz;
        for (k=0;k<3;k++) pdata[i].gcmvel[k] /= EncMass;
    }
    delete[] numingroup;
    t1=MyGetTime()-t1;
    cout<<"Done getting cm in "<<t1<<endl;
}

void MovetoMBFrame(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp){
    cout<<"Moving groups to Most bound particle frame"<<endl;
    Particle *Pval;
    Int_t i,j,k,noffset;
    Double_t x,y,z,vx,vy,vz;
    Double_t t1=MyGetTime();
    //for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,Pval,x,y,z,vx,vy,vz,noffset)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (hp[i].NumberofParticles<=omppropnum)
    {
        noffset=hp[i].noffset;
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - pdata[i].gpos[0];
            y = (*Pval).Y() - pdata[i].gpos[1];
            z = (*Pval).Z() - pdata[i].gpos[2];
            vx = (*Pval).Vx() - pdata[i].gvel[0];
            vy = (*Pval).Vy() - pdata[i].gvel[1];
            vz = (*Pval).Vz() - pdata[i].gvel[2];
            Pval->SetPosition(x,y,z);
            Pval->SetVelocity(vx,vy,vz);
        }
    }
#ifdef USEOPENMP
}
#endif
    for (i=0;i<ngroup;i++) if (hp[i].NumberofParticles>omppropnum)
    {
        noffset=hp[i].noffset;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,Pval,x,y,z,vx,vy,vz)
{
    #pragma omp for 
#endif
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - pdata[i].gpos[0];
            y = (*Pval).Y() - pdata[i].gpos[1];
            z = (*Pval).Z() - pdata[i].gpos[2];
            vx = (*Pval).Vx() - pdata[i].gvel[0];
            vy = (*Pval).Vy() - pdata[i].gvel[1];
            vz = (*Pval).Vz() - pdata[i].gvel[2];
            Pval->SetPosition(x,y,z);
            Pval->SetVelocity(vx,vy,vz);
        }
#ifdef USEOPENMP
}
#endif
    }
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

void MovetoCMFrame(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp){
    cout<<"Moving groups to CM Frame"<<endl;
    Particle *Pval;
    Int_t i,j,k,noffset;
    Double_t x,y,z,vx,vy,vz;
    Double_t t1=MyGetTime();
    //for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,Pval,x,y,z,vx,vy,vz,noffset)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (hp[i].NumberofParticles<=omppropnum)
    {
        noffset=hp[i].noffset;
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - pdata[i].gcm[0];
            y = (*Pval).Y() - pdata[i].gcm[1];
            z = (*Pval).Z() - pdata[i].gcm[2];
            vx = (*Pval).Vx() - pdata[i].gcmvel[0];
            vy = (*Pval).Vy() - pdata[i].gcmvel[1];
            vz = (*Pval).Vz() - pdata[i].gcmvel[2];
            Pval->SetPosition(x,y,z);
            Pval->SetVelocity(vx,vy,vz);
        }
    }
#ifdef USEOPENMP
}
#endif
    for (i=0;i<ngroup;i++) if (hp[i].NumberofParticles>omppropnum)
    {
        noffset=hp[i].noffset;
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,Pval,x,y,z,vx,vy,vz)
{
    #pragma omp for
#endif
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - pdata[i].gcm[0];
            y = (*Pval).Y() - pdata[i].gcm[1];
            z = (*Pval).Z() - pdata[i].gcm[2];
            vx = (*Pval).Vx() - pdata[i].gcmvel[0];
            vy = (*Pval).Vy() - pdata[i].gcmvel[1];
            vz = (*Pval).Vz() - pdata[i].gcmvel[2];
            Pval->SetPosition(x,y,z);
            Pval->SetVelocity(vx,vy,vz);
        }
#ifdef USEOPENMP
}
#endif
    }
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

///move all particles to around particle at deepest point of the potential well
void MovetoPotFrame(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp){
    cout<<"Moving groups to deepest potential particle frame"<<endl;
    GetPotentialEnergy(opt,Part,ngroup,pdata,hp);
    Particle *Pval;
    Int_t i,j,k,noffset,ipot,npot;
    Double_t x,y,z,vx,vy,vz,potval,menc;
    Double_t t1=MyGetTime();
    //use the mean velocity of the max(opt.Nbref,opt.fracbref*numingroup[i]) most bound particles to define the kinetic centre of the group, usually Nbref=10, fracbref=0.05
    //for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,Pval,x,y,z,vx,vy,vz,noffset,menc,ipot,npot,potval)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (hp[i].NumberofParticles<=omppropnum)
    {
        noffset=hp[i].noffset;
        npot=max(opt.Nbref,Int_t(opt.fracbref*hp[i].NumberofParticles));
        potval=MAXVALUE;
        for (j=0;j<hp[i].NumberofParticles;j++) if (Part[noffset+j].GetPotential()<potval) {potval=Part[noffset+j].GetPotential();ipot=j;}
        noffset=hp[i].noffset;
        for (j=0;j<3;j++) {pdata[i].gpos[j]=Part[noffset+ipot].GetPosition(j);pdata[i].gvel[j]=0;}
        menc=0;
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - pdata[i].gpos[0];
            y = (*Pval).Y() - pdata[i].gpos[1];
            z = (*Pval).Z() - pdata[i].gpos[2];
            Pval->SetPosition(x,y,z);
        }
        gsl_heapsort(&Part[noffset],hp[i].NumberofParticles,sizeof(Particle),PotCompare);
        for (j=0;j<npot;j++)
        {
            for (k=0;k<3;k++) pdata[i].gvel[k]+=Part[j+noffset].GetVelocity(k)*Part[j+noffset].GetMass();
            menc+=Part[j+noffset].GetMass();
        }
        pdata[i].gvel=pdata[i].gvel*1.0/menc;
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            vx = (*Pval).Vx() - pdata[i].gvel[0];
            vy = (*Pval).Vy() - pdata[i].gvel[1];
            vz = (*Pval).Vz() - pdata[i].gvel[2];
            Pval->SetVelocity(vx,vy,vz);
        }
        gsl_heapsort(&Part[noffset],hp[i].NumberofParticles,sizeof(Particle),TypeCompare);
    }
#ifdef USEOPENMP
}
#endif
    for (i=0;i<ngroup;i++) if (hp[i].NumberofParticles>omppropnum)
    {
        noffset=hp[i].noffset;
        npot=max(opt.Nbref,Int_t(opt.fracbref*hp[i].NumberofParticles));
        potval=MAXVALUE;
        for (j=0;j<hp[i].NumberofParticles;j++) if (Part[noffset+j].GetPotential()<potval) {potval=Part[noffset+j].GetPotential();ipot=j;}
        for (j=0;j<3;j++) {pdata[i].gpos[j]=Part[noffset+ipot].GetPosition(j);pdata[i].gvel[j]=0;}
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            x = (*Pval).X() - pdata[i].gpos[0];
            y = (*Pval).Y() - pdata[i].gpos[1];
            z = (*Pval).Z() - pdata[i].gpos[2];
            Pval->SetPosition(x,y,z);
        }
        qsort(&Part[noffset],hp[i].NumberofParticles,sizeof(Particle),RadCompare);
        menc=0;
        for (j=0;j<npot;j++) {
            for (k=0;k<3;k++) pdata[i].gvel[k]+=Part[j+noffset].GetVelocity(k)*Part[j+noffset].GetMass();
            menc+=Part[j+noffset].GetMass();
        }
        pdata[i].gvel*=1.0/menc;
/*#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(j,Pval,x,y,z,vx,vy,vz)
{
    #pragma omp for
#endif*/
        for (j=0;j<hp[i].NumberofParticles;j++)
        {
            Pval=&Part[j+noffset];
            vx = (*Pval).Vx() - pdata[i].gvel[0];
            vy = (*Pval).Vy() - pdata[i].gvel[1];
            vz = (*Pval).Vz() - pdata[i].gvel[2];
            Pval->SetVelocity(vx,vy,vz);
        }
/*#ifdef USEOPENMP
}
#endif*/
        qsort(&Part[noffset],hp[i].NumberofParticles,sizeof(Particle),TypeCompare);
    }
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

///Gets bulk properties of groups assuming particles are in group order and the noffsets between groups is given by noffset.
///Assumes that data is sorted according to particle type, then binding energy
///resorts by radius
void GetBulkProp(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp, int ptype)
{
    cout<<"Get Bulk Properties of particles of type "<<ptype<<endl;
    Particle *Pval;
    Int_t i,j,k,noffset,numinvir;
    Double_t ri,rcmv,r2,cmx,cmy,cmz;
    Double_t EncMass,vc,rc,x,y,z,vx,vy,vz;
    Double_t tol=1e-2;
    Int_t ii,icmv;
    Double_t virval=log(opt.virlevel*opt.rhobg);
    Double_t *darray;

    Double_t t1=MyGetTime();

    Int_t *numingroup=new Int_t[ngroup];
    for (i=0;i<ngroup;i++) {
        if (ptype==-1) numingroup[i]=hp[i].NumberofParticles;
        else numingroup[i]=hp[i].NumofType[ptype];
        pdata[i].num=numingroup[i];
    }

    //for small groups loop over groups
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,Pval,ri,rcmv,r2,cmx,cmy,cmz,EncMass,x,y,z,vc,rc,vx,vy,vz,noffset,darray,numinvir,tol)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (numingroup[i]<=omppropnum&&numingroup[i]>0)
    {
        tol=max(1e-2,1.0/sqrt((Double_t)numingroup[i]));
        //set offset correctly
        noffset=hp[i].noffset;
        for (j=1;j<=ptype;j++) noffset+=hp[i].AllNumofType[j-1];
        //sort by radius
        gsl_heapsort(&Part[noffset], numingroup[i], sizeof(Particle), RadCompare);
        pdata[i].gsize=Part[noffset+numingroup[i]-1].Radius();
        //assume in most bound particle frame
        pdata[i].gRmbp=pdata[i].gsize;
        if (pdata[i].gmass==0) for (j=0;j<numingroup[i];j++) pdata[i].gmass+=Part[noffset+j].GetMass();
#ifdef NOMASS
        pdata[i].gmass*=opt.MassValue;
#endif
        //if particles are gas, then also calculate values (means, std, quantiles) characterizing
        //the temperature, metallicty, pressure, etc
        if (numingroup[i]>MINNPROP) {
#ifdef GASON
        if (ptype==GASTYPE) {
            darray=new Double_t[numingroup[i]];
            //first temperature, note that boltzmann stored as log value and output log statistics (equation for temperature is u=3/2*k_b (N_h*Xe+(1-Xe)*N_He)*T
            //all of which is summed up in utoT
            for (j=0;j<numingroup[i];j++) darray[j]=log(Part[j+noffset].GetU())-opt.utoT;
            GetStats(numingroup[i],darray,pdata[i].Temp,opt.nquants,opt.quants);
#ifdef RADIATIVE
            for (j=0;j<numingroup[i];j++) darray[j]=log(Part[j+noffset].GetZmet());
            GetStats(numingroup[i],darray,pdata[i].Zmet,opt.nquants,opt.quants);
#endif
            delete[] darray;
        }
#endif

#ifdef RADIATIVE
        if (ptype==STARTYPE) {
            darray=new Double_t[numingroup[i]];
            //first temperature, note that boltzmann stored as log value and output log statistics
            for (j=0;j<numingroup[i];j++) darray[j]=Part[j+noffset].GetTage();
            GetStats(numingroup[i],darray,pdata[i].Tage,opt.nquants,opt.quants);
            for (j=0;j<numingroup[i];j++) darray[j]=log(Part[j+noffset].GetZmet());
            GetStats(numingroup[i],darray,pdata[i].Zmet,opt.nquants,opt.quants);
            delete[] darray;
        }
#endif
        //determine virial radius
        EncMass=pdata[i].gmass;
        pdata[i].gRvir=pdata[i].gsize;
        numinvir=-1;
        if (ptype==DMTYPE) {
            for (j=numingroup[i]-1;j>=0;j--) {
                Pval=&Part[j+noffset];
                rc=Pval->Radius();
                if (rc>0) 
                    if (-log(4.0*M_PI/3.0)+log(EncMass)-3.0*log(rc)>virval) {pdata[i].gMvir=EncMass;pdata[i].gRvir=rc;numinvir=j+1;break;}
#ifdef NOMASS
                EncMass-=Pval->GetMass();
#else
                EncMass-=opt.MassValue*Pval->GetMass();
#endif
            }
            if (numinvir<0) {pdata[i].gMvir=pdata[i].gmass;pdata[i].gRvir=pdata[i].gsize;numinvir=numingroup[i];}
            else if (numinvir<MINNPROP) numinvir=min((Int_t)numingroup[i],(Int_t)MINNPROP);
        }


        //then calculate all quantities using this region
        pdata[i].gq=pdata[i].gs=1.0;
#ifdef NOMASS
        GetGlobalSpatialMorphology(numingroup[i], &Part[noffset], pdata[i].gq, pdata[i].gs, tol, pdata[i].geigvec,0);
#else
        GetGlobalSpatialMorphology(numingroup[i], &Part[noffset], pdata[i].gq, pdata[i].gs, tol, pdata[i].geigvec,1);
#endif
        pdata[i].gmaxvel=0.;
        EncMass=0;
        pdata[i].gJ[0]=pdata[i].gJ[1]=pdata[i].gJ[2]=0.;
        for (j=0;j<numingroup[i];j++) {
            vc=0;
            Pval=&Part[j+noffset];
            EncMass+=Pval->GetMass();
#ifdef NOMASS
            EncMass*=opt.MassValue;
#endif
            rc=Pval->Radius();
            vx = (*Pval).Vx();
            vy = (*Pval).Vy();
            vz = (*Pval).Vz();
            pdata[i].gJ=pdata[i].gJ+Coordinate(Pval->GetPosition()).Cross(Coordinate(vx,vy,vz))*Pval->GetMass();
            pdata[i].gveldisp(0,0)+=vx*vx*Pval->GetMass();
            pdata[i].gveldisp(1,1)+=vy*vy*Pval->GetMass();
            pdata[i].gveldisp(2,2)+=vz*vz*Pval->GetMass();
            pdata[i].gveldisp(0,1)+=vx*vy*Pval->GetMass();
            pdata[i].gveldisp(0,2)+=vx*vz*Pval->GetMass();
            pdata[i].gveldisp(1,2)+=vy*vz*Pval->GetMass();
            if (rc>0) if (EncMass>0) vc=sqrt(opt.G*EncMass*opt.MassValue/rc);
            //max circ and then vir data
            if (vc>pdata[i].gmaxvel&&EncMass>=MINMASSFRAC*pdata[i].gmass) {pdata[i].gmaxvel=vc;pdata[i].gRmaxvel=rc;pdata[i].gMmaxvel=EncMass;}
        }
        pdata[i].gveldisp(1,0)=pdata[i].gveldisp(0,1);
        pdata[i].gveldisp(2,0)=pdata[i].gveldisp(0,2);
        pdata[i].gveldisp(2,1)=pdata[i].gveldisp(1,2);
        pdata[i].gveldisp=pdata[i].gveldisp*(1.0/EncMass);
        pdata[i].gsigma_v=pow(pdata[i].gveldisp.Det(),1.0/6.0);
#ifdef NOMASS
        pdata[i].gJ=pdata[i].gJ*opt.MassValue;
        pdata[i].gMvir*=opt.MassValue;
        pdata[i].gMmaxvel*=opt.MassValue;
        pdata[i].gmass*=opt.MassValue;
#endif
        }
    }
#ifdef USEOPENMP
}
#endif

    for (i=0;i<ngroup;i++) if (numingroup[i]>omppropnum)
    {
        tol=max(1e-2,1.0/sqrt((Double_t)numingroup[i]));
        //set offset correctly
        noffset=hp[i].noffset;
        for (j=1;j<=ptype;j++) noffset+=hp[i].AllNumofType[j-1];
        qsort(&Part[noffset], numingroup[i], sizeof(Particle), RadCompare);
        pdata[i].gsize=Part[noffset+numingroup[i]-1].Radius();
        pdata[i].gRmbp=pdata[i].gsize;
        if (pdata[i].gmass==0) for (j=0;j<numingroup[i];j++) pdata[i].gmass+=Part[noffset+j].GetMass();
#ifdef NOMASS
        pdata[i].gmass*=opt.MassValue;
#endif
        //if particles are gas, then also calculate values (means, std, quantiles) characterizing
        //the temperature, metallicty, pressure, etc
#ifdef GASON
        if (ptype==GASTYPE) {
            darray=new Double_t[numingroup[i]];
            //first temperature, note that boltzmann stored as log value and output log statistics
            for (j=0;j<numingroup[i];j++) darray[j]=log(Part[j+noffset].GetU())-opt.utoT;
            GetStats(numingroup[i],darray,pdata[i].Temp,opt.nquants,opt.quants);
#ifdef RADIATIVE
            for (j=0;j<numingroup[i];j++) darray[j]=log(Part[j+noffset].GetZmet());
            GetStats(numingroup[i],darray,pdata[i].Zmet,opt.nquants,opt.quants);
#endif
            delete[] darray;
        }
#endif
#ifdef RADIATIVE
        if (ptype==STARTYPE) {
            darray=new Double_t[numingroup[i]];
            //first temperature, note that boltzmann stored as log value and output log statistics
            for (j=0;j<numingroup[i];j++) darray[j]=Part[j+noffset].GetTage();
            GetStats(numingroup[i],darray,pdata[i].Tage,opt.nquants,opt.quants);
            for (j=0;j<numingroup[i];j++) darray[j]=log(Part[j+noffset].GetZmet());
            GetStats(numingroup[i],darray,pdata[i].Zmet,opt.nquants,opt.quants);
            delete[] darray;
        }
#endif
        EncMass=pdata[i].gmass;
        pdata[i].gRvir=0.;
        numinvir=0;
        if (ptype==DMTYPE) {
            for (j=numingroup[i]-1;j>0;j--) {
                Pval=&Part[j+noffset];
                rc=Pval->Radius();
                if (rc>0) 
                    if (-log(4.0*M_PI/3.0)+log(EncMass)-3.0*log(rc)>virval) {pdata[i].gMvir=EncMass;pdata[i].gRvir=rc;numinvir=j+1;break;}
#ifdef NOMASS
                EncMass-=Pval->GetMass();
#else
                EncMass-=opt.MassValue*Pval->GetMass();
#endif
            }
            if (numinvir==0) {pdata[i].gMvir=pdata[i].gmass;pdata[i].gRvir=pdata[i].gsize;numinvir=numingroup[i];}
            else if (numinvir<MINNPROP) numinvir=min((Int_t)numingroup[i],(Int_t)MINNPROP);
        }

        pdata[i].gq=pdata[i].gs=1.0;
#ifdef NOMASS
        GetGlobalSpatialMorphology(numingroup[i], &Part[noffset], pdata[i].gq, pdata[i].gs, tol, pdata[i].geigvec,0);
#else
        GetGlobalSpatialMorphology(numingroup[i], &Part[noffset], pdata[i].gq, pdata[i].gs, tol, pdata[i].geigvec,1);
#endif
        Double_t Jx,Jy,Jz,sxx,sxy,sxz,syy,syz,szz;
        Coordinate J;
        Jx=Jy=Jz=sxx=sxy=sxz=syy=syz=szz=0.;
        EncMass=0;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(j,Pval,x,y,z,vx,vy,vz,J)
{
    #pragma omp for reduction(+:EncMass,Jx,Jy,Jz,sxx,sxy,sxz,syy,syz,szz)
#endif
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            EncMass+=Pval->GetMass();
#ifdef NOMASS
            EncMass*=opt.MassValue;
#endif
            vx = (*Pval).Vx();
            vy = (*Pval).Vy();
            vz = (*Pval).Vz();
            J=Coordinate(Pval->GetPosition()).Cross(Coordinate(vx,vy,vz))*Pval->GetMass();
            Jx+=J[0];Jy+=J[1];Jz+=J[2];
            sxx+=vx*vx*(*Pval).GetMass();
            syy+=vy*vy*(*Pval).GetMass();
            szz+=vz*vz*(*Pval).GetMass();
            sxy+=vx*vy*(*Pval).GetMass();
            sxz+=vx*vz*(*Pval).GetMass();
            syz+=vy*vz*(*Pval).GetMass();
        }
#ifdef USEOPENMP
}
#endif

        pdata[i].gJ[0]=Jx*opt.MassValue;
        pdata[i].gJ[1]=Jy*opt.MassValue;
        pdata[i].gJ[2]=Jz*opt.MassValue;
        pdata[i].gveldisp(0,0)=sxx;
        pdata[i].gveldisp(1,1)=syy;
        pdata[i].gveldisp(2,2)=szz;
        pdata[i].gveldisp(0,1)=pdata[i].gveldisp(1,0)=sxy;
        pdata[i].gveldisp(0,2)=pdata[i].gveldisp(2,0)=sxz;
        pdata[i].gveldisp(1,2)=pdata[i].gveldisp(2,1)=syz;
        pdata[i].gveldisp=pdata[i].gveldisp*(1.0/EncMass);
        pdata[i].gsigma_v=pow(pdata[i].gveldisp.Det(),1.0/6.0);
        EncMass=0.;
        pdata[i].gmaxvel=0.;
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            EncMass+=Pval->GetMass();
#ifdef NOMASS
            EncMass*=opt.MassValue;
#endif
            vc=0;
            rc=Pval->Radius();
            if (rc>0) if (EncMass>0) vc=sqrt(opt.G*EncMass/rc);
            //ensure enough mass in vmax region
            if (vc>pdata[i].gmaxvel && EncMass>MINMASSFRAC*pdata[i].gmass) {pdata[i].gmaxvel=vc;pdata[i].gRmaxvel=rc;pdata[i].gMmaxvel=EncMass;}
        }
#ifdef NOMASS
        pdata[i].gmass*=opt.MassValue;
        pdata[i].gMvir*=opt.MassValue;
        pdata[i].gMmaxvel*=opt.MassValue;
#endif
    }
    delete[] numingroup;
    t1=MyGetTime()-t1;
    cout<<"Done getting properties in "<<t1<<endl;
}

///Get spatial morphology using iterative procedure
void GetGlobalSpatialMorphology(const Int_t nbodies, Particle *p, Double_t& q, Double_t& s, Double_t Error, Matrix& eigenvec, int imflag)
{
    // Calculate the axial ratios q and s.
    Double_t oldq,olds;
    Coordinate e;
    Matrix R(0.),M(0.0),eigenvecp(0.);
    eigenvec=Matrix(0.);
    eigenvec(0,0)=eigenvec(1,1)=eigenvec(2,2)=1.0;
    // Iterative procedure.  See Dubinski and Carlberg (1991).
    int i=0;
    do
    {
        M = Matrix(0.0);
        eigenvecp=Matrix(0.);
        if (imflag==1)CalcMTensorWithMass(M, q, s, nbodies, p);
        else CalcMTensor(M, q, s, nbodies, p);
        e = M.Eigenvalues();
        oldq = q;olds = s;
        q = sqrt(e[1] / e[0]);s = sqrt(e[2] / e[0]);
        eigenvecp=M.Eigenvectors(e);
        eigenvec=eigenvecp*eigenvec;
        RotParticles(nbodies, p, eigenvecp);
        i++;
    } while (fabs((olds-s)/s) >Error || fabs(oldq - q)/q > Error);
    //rotate system back to original coordinate frame
    R=eigenvec.Transpose();
    RotParticles(nbodies, p, R);
}

///calculate the inertia tensor
void CalcITensor(const Int_t n, Particle *p, Double_t &a, Double_t &b, Double_t &c, Matrix& eigenvec, Matrix &I)
{
    Double_t r2,det,Ixx,Iyy,Izz,Ixy,Ixz,Iyz;
    Coordinate e;
    I=Matrix(0.);
    Int_t i;
    Ixx=Iyy=Izz=Ixy=Ixz=Iyz=0.;
#ifdef USEOPENMP
    if (n>omppropnum) {
#pragma omp parallel default(shared) \
private(i,r2)
{
#pragma omp for schedule(dynamic) nowait reduction(+:Ixx,Iyy,Izz,Ixy,Ixz,Iyz)
    for (i = 0; i < n; i++)
    {
        r2=p[i].X()*p[i].X()+p[i].Y()*p[i].Y()+p[i].Z()*p[i].Z();
        Ixx+=p[i].GetMass()*(r2-p[i].X()*p[i].X());
        Iyy+=p[i].GetMass()*(r2-p[i].Y()*p[i].Y());
        Izz+=p[i].GetMass()*(r2-p[i].Z()*p[i].Z());
        Ixy+=p[i].GetMass()*(-p[i].X()*p[i].Y());
        Ixz+=p[i].GetMass()*(-p[i].X()*p[i].Z());
        Iyz+=p[i].GetMass()*(-p[i].Y()*p[i].Z());
    }
}
    I(0,0)=Ixx;I(1,1)=Iyy;I(2,2)=Izz;
    I(0,1)=I(1,0)=Ixy;
    I(0,2)=I(2,0)=Ixz;
    I(1,2)=I(2,1)=Iyz;
    }
    else {
#endif
    for (i = 0; i < n; i++)
    {
        r2=p[i].X()*p[i].X()+p[i].Y()*p[i].Y()+p[i].Z()*p[i].Z();
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                I(j, k) += p[i].GetMass()*((j==k)*r2-p[i].GetPosition(j)*p[i].GetPosition(k));
            }
        }
    }
#ifdef USEOPENMP
    }
#endif
    det=I.Det();
    for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) I(j, k) /= det;
    e = I.Eigenvalues();
    a=e[0];b=e[1];c=e[2];
    eigenvec=I.Eigenvectors(e);
}

///calculate the weighted reduced inertia tensor assuming particles are the same mass
void CalcMTensor(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p)
{
    Int_t i;
    int j,k;
    Double_t a2,Mxx,Myy,Mzz,Mxy,Mxz,Myz;
    Mxx=Myy=Mzz=Mxy=Mxz=Myz=0.;
#ifdef USEOPENMP
    if (n>omppropnum) {
#pragma omp parallel default(shared) \
private(i,a2)
{
#pragma omp for schedule(dynamic) nowait reduction(+:Mxx,Myy,Mzz,Mxy,Mxz,Myz)
    for (i = 0; i < n; i++)
    {
        a2 = (p[i].X()*p[i].X()+p[i].Y()*p[i].Y()/q/q+p[i].Z()*p[i].Z()/s/s);
        if (a2>0) {
        a2=1.0/a2;
        Mxx+=p[i].X()*p[i].X()*a2;
        Myy+=p[i].Y()*p[i].Y()*a2;
        Mzz+=p[i].Z()*p[i].Z()*a2;
        Mxy+=p[i].X()*p[i].Y()*a2;
        Mxz+=p[i].X()*p[i].Z()*a2;
        Myz+=p[i].Y()*p[i].Z()*a2;
        }
    }
}
    M(0,0)=Mxx;M(1,1)=Myy;M(2,2)=Mzz;
    M(0,1)=M(1,0)=Mxy;
    M(0,2)=M(2,0)=Mxz;
    M(1,2)=M(2,1)=Myz;
    }
    else {
#endif
    for (i = 0; i < n; i++)
    {
        a2 = (p[i].X()*p[i].X()+p[i].Y()*p[i].Y()/q/q+p[i].Z()*p[i].Z()/s/s);
        if(a2>0){
            a2=1.0/a2;
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                M(j, k) += p[i].GetPosition(j)*p[i].GetPosition(k)*a2;
        }
    }
#ifdef USEOPENMP
    }
#endif
}

///calculate the weighted reduced inertia tensor
void CalcMTensorWithMass(Matrix& M, const Double_t q, const Double_t s, const Int_t n, Particle *p)
{
    Int_t i;
    int j,k;
    Double_t a2,Mxx,Myy,Mzz,Mxy,Mxz,Myz;
    Mxx=Myy=Mzz=Mxy=Mxz=Myz=0.;
#ifdef USEOPENMP
    if (n>omppropnum) {
#pragma omp parallel default(shared) \
private(i,a2)
{
#pragma omp for schedule(dynamic) nowait reduction(+:Mxx,Myy,Mzz,Mxy,Mxz,Myz)
    for (i = 0; i < n; i++)
    {
        a2 = (p[i].X()*p[i].X()+p[i].Y()*p[i].Y()/q/q+p[i].Z()*p[i].Z()/s/s);
        if (a2>0) {
        a2=p[i].GetMass()/a2;
        Mxx+=p[i].X()*p[i].X()*a2;
        Myy+=p[i].Y()*p[i].Y()*a2;
        Mzz+=p[i].Z()*p[i].Z()*a2;
        Mxy+=p[i].X()*p[i].Y()*a2;
        Mxz+=p[i].X()*p[i].Z()*a2;
        Myz+=p[i].Y()*p[i].Z()*a2;
        }
    }
}
    M(0,0)=Mxx;M(1,1)=Myy;M(2,2)=Mzz;
    M(0,1)=M(1,0)=Mxy;
    M(0,2)=M(2,0)=Mxz;
    M(1,2)=M(2,1)=Myz;
    }
    else {
#endif
    for (i = 0; i < n; i++)
    {
        a2 = (p[i].X()*p[i].X()+p[i].Y()*p[i].Y()/q/q+p[i].Z()*p[i].Z()/s/s);
        if(a2>0){
            a2=p[i].GetMass()/a2;
            for (j = 0; j < 3; j++)
                for (k = 0; k < 3; k++)
                M(j, k) += p[i].GetPosition(j)*p[i].GetPosition(k)*a2;
        }
    }
#ifdef USEOPENMP
    }
#endif
}

///rotate particles
void RotParticles(const Int_t n, Particle *p, Matrix &R)
{
    Int_t i;
    int j;
    Coordinate temp;
#ifdef USEOPENMP
    if (n>omppropnum) {
#pragma omp parallel default(shared) \
private(i,j,temp)
{
#pragma omp for schedule(dynamic) nowait
    for (i=0; i<n; i++)
    {
        temp[0]=temp[1]=temp[2]=0.;
        for (j=0; j<3; j++)
        {
            temp[0] += R(0, j) * (p[i].GetPosition(j));
            temp[1] += R(1, j) * (p[i].GetPosition(j));
            temp[2] += R(2, j) * (p[i].GetPosition(j));
        }
        p[i].SetPosition(temp[0], temp[1], temp[2]);
    }
}
    }
    else {
#endif
    for (i=0; i<n; i++)
    {
        temp[0]=temp[1]=temp[2]=0.;
        for (j=0; j<3; j++)
        {
            temp[0] += R(0, j) * (p[i].GetPosition(j));
            temp[1] += R(1, j) * (p[i].GetPosition(j));
            temp[2] += R(2, j) * (p[i].GetPosition(j));
        }
        p[i].SetPosition(temp[0], temp[1], temp[2]);
    }
#ifdef USEOPENMP
    }
#endif
}


/*!
    Get various radial profiles
    Assumes that data is sorted type, then radius and that total energy is stored in GetPotential()
    one can sort via energy, determine a cut off point, resort by radius within that energy cutoff and
    then call this routine to get the radial profile of bound particles. One just needs to update
    the appropriate \ref NumberofParticles of \ref NumofType in \ref HaloParticleData
    For DM:
        Menc
        q,s
        eigvec
        velocity dispersion tensor
    For gas, same but also
        mass weighted temperature (based on self energy where U=(3/2)*(k_B)*T)
        density (based on sph density)
        pressure (based on sph density, \f$ P=A\rho^{\gamma}\f$ where A is the entropy
        mass weighted metallicity
    From these quantities, one can calculate a variety of quantities such star formation rate,
    as a mock X-ray radial profile which can be projected to a surface using the eigenvectors of the mass
    distribution
    For stars, same as dm but also
        Mass weighted age
        mass weighted metallicity
    Its possible that will add more. Also will add other specific calculations valid for gas or stars
*/
void GetProfiles(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp, int ptype){
    cout<<"Getting Profiles "<<endl;
    Particle *Pval;
    Int_t i,j,k,n,noffset;
    Int_t nbins,npartbin,ninbin,ndatablocks,*idatablocks;
    Double_t norm,tol;
    Double_t t1;
    Double_t vx,vy,vz,sxx,sxy,sxz,syy,syz,szz;
    Double_t logTmean,Temweight,Tslweight;
    //set the data blocks to be analyzed
    if (ptype==0) {
        ndatablocks=5;
    }
    t1=MyGetTime();

    Int_t *numingroup=new Int_t[ngroup];
    for (i=0;i<ngroup;i++) {
        if (ptype==-1) numingroup[i]=hp[i].NumberofParticles;
        else numingroup[i]=hp[i].NumofType[ptype];
    }

#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,n,Pval,noffset,nbins,npartbin,ninbin,norm,vx,vy,vz,sxx,sxy,sxz,syy,syz,szz,logTmean,Temweight,Tslweight,tol)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif

    for (i=0;i<ngroup;i++) {
    if (numingroup[i]>MINNPROFILE)
    {
        noffset=hp[i].noffset;
        if (ptype>0) for (j=1;j<=ptype;j++) noffset+=hp[i].AllNumofType[j-1];
        //have bin widths vary and contain approximately the same number of particles
        nbins=GetNBins(numingroup[i]);
        npartbin=ceil(numingroup[i]/(Double_t)nbins);

        //pdata[i].radprofile.nbins=nbins;
        pdata[i].radprofile.SetBins(nbins);

        //since sorted radially, don't have to check radius but might have to check energy
        //issue with this choice of radial bins is that different particle types will have different
        //radial bins. Just means one has to interpolate to get gas fraction as a function of radius
        for (k=0;k<nbins;k++) {
            pdata[i].radprofile.Menc[k]=0;
            pdata[i].radprofile.Jenc[k][0]=pdata[i].radprofile.Jenc[k][1]=pdata[i].radprofile.Jenc[k][2]=0.;
            pdata[i].radprofile.radval[k]=0;
            pdata[i].radprofile.qenc[k]=pdata[i].radprofile.senc[k]=1.0;
            if (ptype==GASTYPE) {
                pdata[i].radprofile.Tempenc[k]=0.;
#ifdef RADIATIVE
                pdata[i].radprofile.Zmetenc[k]=0.;
#endif
            }
#ifdef RADIATIVE
            if (ptype==STARTYPE) {
                pdata[i].radprofile.Zmetenc[k]=0.;
                pdata[i].radprofile.Tageenc[k]=0.;
            }
#endif
            if (k!=nbins-1) n=(k+1)*npartbin;
            else n=numingroup[i];
            ninbin=n-k*npartbin;
            for (j=k*npartbin;j<n;j++)
            {
                Pval=&Part[j+noffset];
                pdata[i].radprofile.numinbin[k]+=1.0;
                pdata[i].radprofile.radval[k]+=Pval->Radius()*Pval->GetMass();
                pdata[i].radprofile.Menc[k]+=Pval->GetMass();
                pdata[i].radprofile.Jenc[k][0]+=Pval->GetMass()*(Pval->Y()*Pval->Vz()-Pval->Vy()*Pval->Z());
                pdata[i].radprofile.Jenc[k][1]+=Pval->GetMass()*(-Pval->X()*Pval->Vz()+Pval->Vx()*Pval->Z());
                pdata[i].radprofile.Jenc[k][2]+=Pval->GetMass()*(Pval->X()*Pval->Vy()-Pval->Vx()*Pval->Y());
                for (int kk=0;kk<3;kk++) pdata[i].radprofile.velenc[k][kk]+=Pval->GetMass()*Pval->GetVelocity(kk);
            }
            pdata[i].radprofile.radval[k]/=pdata[i].radprofile.Menc[k];
            for (int kk=0;kk<3;kk++) pdata[i].radprofile.velenc[k][kk]/=pdata[i].radprofile.Menc[k];
            sxx=syy=szz=sxy=sxz=syz=0.;
            for (j=k*npartbin;j<n;j++)
            {
                Pval=&Part[j+noffset];
                vx = (*Pval).Vx();
                vy = (*Pval).Vy();
                vz = (*Pval).Vz();
                sxx+=vx*vx*(*Pval).GetMass();
                syy+=vy*vy*(*Pval).GetMass();
                szz+=vz*vz*(*Pval).GetMass();
                sxy+=vx*vy*(*Pval).GetMass();
                sxz+=vx*vz*(*Pval).GetMass();
                syz+=vy*vz*(*Pval).GetMass();
            }
            pdata[i].radprofile.sigmavenc[k](0,0)=sxx;
            pdata[i].radprofile.sigmavenc[k](1,1)=syy;
            pdata[i].radprofile.sigmavenc[k](2,2)=szz;
            pdata[i].radprofile.sigmavenc[k](0,1)=pdata[i].radprofile.sigmavenc[k](1,0)=sxy;
            pdata[i].radprofile.sigmavenc[k](0,2)=pdata[i].radprofile.sigmavenc[k](2,0)=sxz;
            pdata[i].radprofile.sigmavenc[k](1,2)=pdata[i].radprofile.sigmavenc[k](2,1)=syz;
            pdata[i].radprofile.sigmavenc[k]=pdata[i].radprofile.sigmavenc[k]*(1.0/pdata[i].radprofile.Menc[k]);
            pdata[i].radprofile.qenc[k]=pdata[i].radprofile.senc[k]=1.0;
            tol=max(1e-2,1.0/sqrt((Double_t)npartbin));
            GetGlobalSpatialMorphology(ninbin, &Part[noffset+k*npartbin], pdata[i].radprofile.qenc[k], pdata[i].radprofile.senc[k], tol, pdata[i].radprofile.eigvecenc[k]);
#ifdef GASON
            if (ptype==GASTYPE) {
            logTmean=Temweight=Tslweight=0;
            for (j=k*npartbin;j<n;j++) {
                Pval=&Part[j+noffset];
                pdata[i].radprofile.Den[k]+=Pval->GetDensity()*Pval->GetMass();
                pdata[i].radprofile.Tempenc[k]+=exp(log(Pval->GetU())-opt.utoT)*Pval->GetMass();
                logTmean+=(log(Pval->GetU())-opt.utoT)*Pval->GetMass();
                pdata[i].radprofile.Tdispenc[k]+=(log(Pval->GetU())-opt.utoT)*(log(Pval->GetU())-opt.utoT)*Pval->GetMass();
                pdata[i].radprofile.Temenc[k]+=exp(log(Pval->GetU())-opt.utoT)*Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                Temweight+=Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                if (log(Pval->GetU())-opt.utoT>opt.Tempthreshold) {
                pdata[i].radprofile.Tslenc[k]+=pow(exp(log(Pval->GetU())-opt.utoT),0.25)*Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                Tslweight+=pow(exp(log(Pval->GetU())-opt.utoT),-0.75)*Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                }
#ifdef RADIATIVE
                pdata[i].radprofile.Zmetenc[k]+=Pval->GetZmet()*Pval->GetMass();
#endif
            }
            pdata[i].radprofile.Den[k]/=pdata[i].radprofile.Menc[k];
            pdata[i].radprofile.Tempenc[k]/=pdata[i].radprofile.Menc[k];
            pdata[i].radprofile.Tdispenc[k]/=pdata[i].radprofile.Menc[k];
            logTmean/=pdata[i].radprofile.Menc[k];
            pdata[i].radprofile.Tdispenc[k]-=(logTmean*logTmean);
            pdata[i].radprofile.Tdispenc[k]=sqrt(pdata[i].radprofile.Tdispenc[k]);
            pdata[i].radprofile.Temenc[k]/=Temweight;
            pdata[i].radprofile.Tslenc[k]/=Tslweight;
#ifdef RADIATIVE
            pdata[i].radprofile.Zmetenc[k]/=pdata[i].radprofile.Menc[k];
#endif
            }
#endif
#ifdef RADIATIVE
            if (ptype==STARTYPE) {
            for (j=k*npartbin;j<n;j++) {
                Pval=&Part[j+noffset];
                pdata[i].radprofile.Tageenc[k]+=Pval->GetTage()*Pval->GetMass();
                pdata[i].radprofile.Zmetenc[k]+=Pval->GetZmet()*Pval->GetMass();
            }
            pdata[i].radprofile.Tageenc[k]/=pdata[i].radprofile.Menc[k];
            pdata[i].radprofile.Zmetenc[k]/=pdata[i].radprofile.Menc[k];
            }
#endif
        }
    }
    else pdata[i].radprofile.nbins=0;
    }
#ifdef USEOPENMP
}
#endif
    delete[] numingroup;
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

void GetCylFrame(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp){
    cout<<"Getting Cylindrical Profiles "<<endl;
    Particle *Pval;
    Int_t i,j,k,n,noffset,noffsettype;
    Int_t nbins,npartbin,ninbin,ndatablocks,*idatablocks;
    Double_t norm;
    Double_t t1;
    Double_t Rval,vx,vy,vz,sxx,sxy,sxz,syy,syz,szz;
    Coordinate newx;
    Matrix Rot;
    t1=MyGetTime();
    Int_t *numingroup=new Int_t[ngroup];
    for (i=0;i<ngroup;i++) {
        numingroup[i]=hp[i].NumberofParticles;
    }
    //prior to this point, particles where sorted by radius. But rotate into angular momentum frame and sort by cylindrical radius
    //rotate into coordinate frame
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,n,Pval,noffset,Rot,norm,newx)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) if (numingroup[i]>0)
    {
        noffset=hp[i].noffset;
        pdata[i].cylprofile.zref=pdata[i].gJ*(1.0/pdata[i].gJ.Length());
        //if not aligned to z-axis already, rotate
        if (!(pdata[i].cylprofile.zref[0]==0&&pdata[i].cylprofile.zref[1]==0&&pdata[i].cylprofile.zref[0]==1)) {
        Rot(2,0)=pdata[i].cylprofile.zref[0];
        Rot(2,1)=pdata[i].cylprofile.zref[1];
        Rot(2,2)=pdata[i].cylprofile.zref[2];
        norm=1.0/sqrt(pdata[i].cylprofile.zref[1]*pdata[i].cylprofile.zref[1]+pdata[i].cylprofile.zref[0]*pdata[i].cylprofile.zref[0]);
        Rot(0,0)=pdata[i].cylprofile.zref[1]*norm;
        Rot(0,1)=-pdata[i].cylprofile.zref[0]*norm;
        Rot(0,2)=0;
        Rot(1,0)=pdata[i].cylprofile.zref[0]*pdata[i].cylprofile.zref[2];
        Rot(1,1)=pdata[i].cylprofile.zref[1]*pdata[i].cylprofile.zref[2];
        Rot(1,2)=-pdata[i].cylprofile.zref[0]*pdata[i].cylprofile.zref[0]-pdata[i].cylprofile.zref[1]*pdata[i].cylprofile.zref[1];
        norm=1.0/sqrt(Rot(1,0)*Rot(1,0)+Rot(1,1)*Rot(1,1)+Rot(1,2)*Rot(1,2));
        Rot(1,0)*=norm;
        Rot(1,1)*=norm;
        Rot(1,2)*=norm;
        for (j=0;j<numingroup[i];j++) {
            Pval=&Part[j+noffset];
            for (k=0;k<3;k++) newx[k]=Pval->GetPosition(k);
            newx=Rot*newx;
            for (k=0;k<3;k++) Pval->SetPosition(k,newx[k]);
            for (k=0;k<3;k++) newx[k]=Pval->GetVelocity(k);
            newx=Rot*newx;
            for (k=0;k<3;k++) Pval->SetVelocity(k,newx[k]);
        }
        }
        //qsort(&Part[noffset], numingroup[i], sizeof(Particle), CylRadCompare);
    }
#ifdef USEOPENMP
}
#endif

/*
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
*/

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
                qsort(&Part[hp[i].noffset+noffsettype], hp[i].NumofType[j], sizeof(Particle), CylRadCompare);
            noffsettype+=hp[i].NumofType[j];
        }
    }
#ifdef USEOPENMP
}
#endif
    delete[] numingroup;
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

void GetCylProfiles(Options &opt, Particle *Part, Int_t ngroup, PropData *pdata, HaloParticleData *hp, int ptype){
    cout<<"Getting Cylindrical Profiles "<<endl;
    Particle *Pval;
    Int_t i,j,k,n,noffset;
    Int_t nbins,npartbin,ninbin,ndatablocks,*idatablocks;
    Double_t norm,tol;
    Double_t t1;
    Double_t Rval,vx,vy,vz,sxx,sxy,sxz,syy,syz,szz;
    Double_t logTmean,Temweight,Tslweight;
    Coordinate newx;
    Matrix Rot;
    //set the data blocks to be analyzed
    if (ptype==0) {
        ndatablocks=5;
    }
    t1=MyGetTime();
    Int_t *numingroup=new Int_t[ngroup];
    for (i=0;i<ngroup;i++) {
        if (ptype==-1) numingroup[i]=hp[i].NumberofParticles;
        else numingroup[i]=hp[i].NumofType[ptype];
    }

#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i,j,k,n,Pval,noffset,nbins,npartbin,ninbin,norm,vx,vy,vz,sxx,sxy,sxz,syy,syz,szz,Rval,logTmean,Temweight,Tslweight,tol)
{
    #pragma omp for schedule(dynamic,1) nowait
#endif
    for (i=0;i<ngroup;i++) {
    if (numingroup[i]>MINNPROFILE) 
    {
        noffset=hp[i].noffset;
        if (ptype>0) for (j=1;j<=ptype;j++) noffset+=hp[i].AllNumofType[j-1];
        nbins=GetNBins(numingroup[i]);
        //have bin widths vary and contain approximately the same number of particles
        npartbin=round(numingroup[i]/(Double_t)nbins);
        //pdata[i].cylprofile.nbins=nbins;
        pdata[i].cylprofile.SetBins(nbins);
        //since sorted radially, don't have to check radius but might have to check energy
        //issue with this choice of radial bins is that different particle types will have different
        //radial bins. Just means one has to interpolate to get gas fraction as a function of radius
        for (k=0;k<nbins;k++) {
            pdata[i].cylprofile.Menc[k]=0;
	        pdata[i].cylprofile.Rmean[k]=pdata[i].cylprofile.Rdisp[k]=pdata[i].cylprofile.zmean[k]=pdata[i].cylprofile.zdisp[k]=0.;
            pdata[i].cylprofile.Jenc[k][0]=pdata[i].cylprofile.Jenc[k][1]=pdata[i].cylprofile.Jenc[k][2]=0.;
            pdata[i].cylprofile.radval[k]=0;
            pdata[i].cylprofile.qenc[k]=pdata[i].cylprofile.senc[k]=1.0;
            if (ptype==GASTYPE) {
                pdata[i].cylprofile.Tempenc[k]=0;
#ifdef RADIATIVE
                pdata[i].cylprofile.Zmetenc[i]=0.;
#endif
            }
#ifdef RADIATIVE
            if (ptype==STARTYPE) {
                pdata[i].cylprofile.Zmetenc[i]=0.;
                pdata[i].cylprofile.Tageenc[k]=0.;
            }
#endif
            n=min((Int_t)((k+1)*npartbin),(Int_t)numingroup[i]);
            ninbin=n-k*npartbin;
            for (j=k*npartbin;j<n;j++)
            {
                Pval=&Part[j+noffset];
                pdata[i].cylprofile.numinbin[k]+=1.0;
                pdata[i].cylprofile.radval[k]+=Pval->Radius()*Pval->GetMass();
                pdata[i].cylprofile.Menc[k]+=Pval->GetMass();
                pdata[i].cylprofile.Jenc[k][0]+=Pval->GetMass()*(Pval->Y()*Pval->Vz()-Pval->Vy()*Pval->Z());
                pdata[i].cylprofile.Jenc[k][1]+=Pval->GetMass()*(-Pval->X()*Pval->Vz()+Pval->Vx()*Pval->Z());
                pdata[i].cylprofile.Jenc[k][2]+=Pval->GetMass()*(Pval->X()*Pval->Vy()-Pval->Vx()*Pval->Y());
                for (int kk=0;kk<3;kk++) pdata[i].cylprofile.velenc[k][kk]+=Pval->GetMass()*Pval->GetVelocity(kk);
                pdata[i].cylprofile.zmean[k]+=Pval->GetMass()*Pval->Z();
                pdata[i].cylprofile.zdisp[k]+=Pval->GetMass()*(Pval->Z()*Pval->Z());
                Rval=Pval->CylRadius();
                pdata[i].cylprofile.Rmean[k]+=Pval->GetMass()*Rval;
            }
            pdata[i].cylprofile.radval[k]/=pdata[i].cylprofile.Menc[k];
            for (int kk=0;kk<3;kk++) pdata[i].cylprofile.velenc[k][kk]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.zmean[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.Rmean[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.zdisp[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.zdisp[k]-=(pdata[i].cylprofile.zmean[k])*(pdata[i].cylprofile.zmean[k]);
            pdata[i].cylprofile.zdisp[k]=sqrt(pdata[i].cylprofile.zdisp[k]);
            for (j=k*npartbin;j<n;j++)
            {
                Pval=&Part[j+noffset];
                Rval=Pval->CylRadius();
                pdata[i].cylprofile.Rdisp[k]+=Pval->GetMass()*pow((Rval-pdata[i].cylprofile.Rmean[k]),2.0);
            }
            pdata[i].cylprofile.Rdisp[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.Rdisp[k]=sqrt(pdata[i].cylprofile.Rdisp[k]);

            sxx=syy=szz=sxy=sxz=syz=0.;
            for (j=k*npartbin;j<n;j++)
            {
                Pval=&Part[j+noffset];
                vx = (*Pval).Vx();
                vy = (*Pval).Vy();
                vz = (*Pval).Vz();
                sxx+=vx*vx*(*Pval).GetMass();
                syy+=vy*vy*(*Pval).GetMass();
                szz+=vz*vz*(*Pval).GetMass();
                sxy+=vx*vy*(*Pval).GetMass();
                sxz+=vx*vz*(*Pval).GetMass();
                syz+=vy*vz*(*Pval).GetMass();
            }
            pdata[i].cylprofile.sigmavenc[k](0,0)=sxx;
            pdata[i].cylprofile.sigmavenc[k](1,1)=syy;
            pdata[i].cylprofile.sigmavenc[k](2,2)=szz;
            pdata[i].cylprofile.sigmavenc[k](0,1)=pdata[i].cylprofile.sigmavenc[k](1,0)=sxy;
            pdata[i].cylprofile.sigmavenc[k](0,2)=pdata[i].cylprofile.sigmavenc[k](2,0)=sxz;
            pdata[i].cylprofile.sigmavenc[k](1,2)=pdata[i].cylprofile.sigmavenc[k](2,1)=syz;
            pdata[i].cylprofile.sigmavenc[k]=pdata[i].cylprofile.sigmavenc[k]*(1.0/pdata[i].cylprofile.Menc[k]);
            pdata[i].cylprofile.qenc[k]=pdata[i].cylprofile.senc[k]=1.0;
            tol=max(1e-2,1.0/sqrt((Double_t)npartbin));
            GetGlobalSpatialMorphology(ninbin, &Part[noffset+k*npartbin], pdata[i].cylprofile.qenc[k], pdata[i].cylprofile.senc[k], tol, pdata[i].cylprofile.eigvecenc[k]);
#ifdef GASON
            if (ptype==GASTYPE) {
            logTmean=Temweight=Tslweight=0;
            for (j=k*npartbin;j<n;j++) {
                Pval=&Part[j+noffset];
                pdata[i].cylprofile.Den[k]+=Pval->GetDensity()*Pval->GetMass();
                pdata[i].cylprofile.Tempenc[k]+=exp(log(Pval->GetU())-opt.utoT)*Pval->GetMass();
                logTmean+=(log(Pval->GetU())-opt.utoT)*Pval->GetMass();
                pdata[i].cylprofile.Tdispenc[k]+=(log(Pval->GetU())-opt.utoT)*(log(Pval->GetU())-opt.utoT)*Pval->GetMass();
                pdata[i].cylprofile.Temenc[k]+=exp(log(Pval->GetU())-opt.utoT)*Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                Temweight+=Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                if (log(Pval->GetU())-opt.utoT>opt.Tempthreshold) {
                pdata[i].cylprofile.Tslenc[k]+=pow(exp(log(Pval->GetU())-opt.utoT),0.25)*Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                Tslweight+=pow(exp(log(Pval->GetU())-opt.utoT),-0.75)*Pval->GetDensity()*Pval->GetDensity()*Pval->GetMass();
                }
#ifdef RADIATIVE
                pdata[i].cylprofile.Zmetenc[k]+=Pval->GetZmet()*Pval->GetMass();
#endif
            }
            pdata[i].cylprofile.Den[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.Tempenc[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.Tdispenc[k]/=pdata[i].cylprofile.Menc[k];
            logTmean/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.Tdispenc[k]-=(logTmean*logTmean);
            pdata[i].cylprofile.Tdispenc[k]=sqrt(pdata[i].cylprofile.Tdispenc[k]);
            pdata[i].cylprofile.Temenc[k]/=Temweight;
            pdata[i].cylprofile.Tslenc[k]/=Tslweight;
#ifdef RADIATIVE
            pdata[i].cylprofile.Zmetenc[k]/=pdata[i].cylprofile.Menc[k];
#endif
            }
#endif
#ifdef RADIATIVE
            if (ptype==STARTYPE) {
            for (j=k*npartbin;j<n;j++) {
                Pval=&Part[j+noffset];
                pdata[i].cylprofile.Tageenc[k]+=Pval->GetTage()*Pval->GetMass();
                pdata[i].cylprofile.Zmetenc[k]+=Pval->GetZmet()*Pval->GetMass();
            }
            pdata[i].cylprofile.Tageenc[k]/=pdata[i].cylprofile.Menc[k];
            pdata[i].cylprofile.Zmetenc[k]/=pdata[i].cylprofile.Menc[k];
            }
#endif
        }
    }
    else pdata[i].cylprofile.nbins=0;
    }
#ifdef USEOPENMP
}
#endif
    delete[] numingroup;
    t1=MyGetTime()-t1;
    cout<<"Done in "<<t1<<endl;
}

int CylRadCompare (const void *a, const void *b){
    Double_t aa = ((Particle*)a)->CylRadius();
    Double_t bb = ((Particle*)b)->CylRadius();
    if (aa > bb) return 1;
    else if (aa < bb) return -1;
    else return 0;
}

inline int GetNBins(Int_t n, int ibintype){
    //sturge's formula modified slightly
    if (ibintype==0) return min((Double_t)NRAD,(Double_t)floor(2.0*(log((Double_t)n)/log(2.0)+1)));
    //just number to inverse of dimension
    else if (ibintype==1) return min((Double_t)NRAD,(Double_t)(pow((Double_t)n,(Double_t)(1./3.))));
}
