/*! \file System.cxx
 *  \brief This file implements functions of \ref NBody::System class

*/


#include <System.h>

using namespace std;
using namespace Math;

namespace NBody
{
    /*-----------------------
        System functions
      -----------------------*/

    // Constructors
    System::System(Int_t num, Double_t t, Double_t *p, Int_t ng, Int_t ns)
    {
        time = t;
        numparts = num;
        ngas=ng;
        nstar=ns;
        ndark=numparts-ngas-nstar;
        particle = new Particle[ndark];
        if (p!=NULL) period=p;
        else for (int i=0.;i<3;i++)period[i]=0.;
        if (ng>0) gparticle = new GasParticle[ng];
        if (ns>0) sparticle = new StarParticle[ns];
    }

    System::System(Int_t num, Particle *P, Double_t t, Double_t *p, Int_t ng, Int_t ns, GasParticle *gp, StarParticle *sp)
    {
        time = t;
        numparts = num;
        ngas=ng;
        nstar=ns;
        ndark=numparts-ngas-nstar;
        particle = P;
        if (p!=NULL) period=p;
        else for (int i=0.;i<3;i++)period[i]=0.;
        if (ng>0) gparticle = gp;
        if (ns>0) sparticle = sp;
    }

    // Copy constructor
    System::System(const System &s)
    {
        numparts = s.numparts;
        time = s.time;
        period = s.period;
        ndark = s.ndark;
        ngas = s.ngas;
        nstar = s.nstar;
        particle = new Particle[ndark];
        for (Int_t i = 0; i < ndark; i++) particle[i] = s.particle[i];
        if (ngas>0) {
            gparticle = new GasParticle[ngas];
            for (Int_t i = 0; i < ngas; i++) gparticle[i] = s.gparticle[i];
        }
        if (nstar>0) {
            sparticle = new StarParticle[nstar];
            for (Int_t i = 0; i < nstar; i++) sparticle[i] = s.sparticle[i];
        }
    }

    //    OTHER FUNCTIONS

    //Properties of the system

    Double_t System::TotalMass() const
    {
      Double_t tm=0;
      for (Int_t i=0;i<numparts;i++)
        tm+=particle[i].GetMass();
      return tm;
    }

    Double_t System::MaxLength() const
    {
      Double_t temp, l=0;
      for (Int_t i=0;i<numparts;i++)
      {
        if (fabs(Part(i).X())>fabs(Part(i).Y()))
        {
          if (fabs(Part(i).X())>fabs(Part(i).Z()))
            temp=fabs(Part(i).X());
          else
            temp=fabs(Part(i).Z());
        }
        else if (fabs(Part(i).Y())>fabs(Part(i).Z()))
        {
          temp=fabs(Part(i).Y());
        }
        else
        {
          temp=fabs(Part(i).Z());
        }

        if(l<temp)
          l=temp;
      }
      return 2.0*l;
    }

    //returns the radius of most distant particle from the CM
    Double_t System::MaxRadius(bool cmframe) const
    {
        Double_t temp, r;
        if (!cmframe) {
            System S(*this);
            S.AdjustForCM();
            r=S[0].Radius();
            for (Int_t i=1;i<numparts;i++)
            {
                temp =S[i].Radius();
                if (temp>r)r=temp;
            }
        }
        else {
            r=particle[0].Radius();
            for (Int_t i=1;i<numparts;i++)
            {
                temp =particle[i].Radius();
                if (temp>r)r=temp;
            }
        }
        return r;
    }

    Double_t System::MinRadius(bool cmframe) const
    {
        Double_t temp, r;
        if (!cmframe) {
            System S(*this);
            S.AdjustForCM();
            r=S[0].Radius();
            for (Int_t i=1;i<numparts;i++)
            {
                temp =S[i].Radius();
                if (temp<r)r=temp;
            }
        }
        else {
            r=particle[0].Radius();
            for (Int_t i=1;i<numparts;i++)
            {
                temp =particle[i].Radius();
                if (temp<r)r=temp;
            }
        }
        return r;
    }

    Double_t *System::RadialLimits(bool cmframe) const
    {
        Double_t temp,*r=new Double_t[3];
        Double_t mtot=0.;
        Double_t averad=0.;
        if (!cmframe) {
            System S(*this);
            S.AdjustForCM();
            r[0]=S[0].Radius();
            r[1]=r[0]*S[0].GetMass();
            r[2]=r[0];
            averad=r[1];
            mtot=S[0].GetMass();

#ifndef USEOPENMP
            for (Int_t i=1;i<numparts;i++)
            {
                temp =S[i].Radius();
                if (temp<r[0])r[0]=temp;
                if (temp>r[2])r[2]=temp;
                r[1]+=temp*S[i].GetMass();
                mtot+=S[i].GetMass();
            }
#else
            Int_t i;
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:averad,mtot)
            for (i=1;i<numparts;i++)
            {
                averad+=S[i].Radius()*S[i].GetMass();
                mtot+=S[i].GetMass();
            }
    }
            r[1]=averad;
            int nthreads;
    #pragma omp parallel 
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
            Double_t *mina=new Double_t[nthreads];
            Double_t *maxa=new Double_t[nthreads];
            for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=r[0];
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic,1) nowait
            for (i=1;i<numparts;i++)
            {
                temp =S[i].Radius();
                if (temp<mina[omp_get_thread_num()])mina[omp_get_thread_num()]=temp;
                if (temp>maxa[omp_get_thread_num()])maxa[omp_get_thread_num()]=temp;
            }
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < r[0]) r[0] = mina[i];
            if (maxa[i] > r[2]) r[2] = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
#endif
        }
        else {
            r[0]=particle[0].Radius();
            r[1]=r[0]*particle[0].GetMass();
            r[2]=r[0];
            averad=r[1];
            mtot=particle[0].GetMass();
#ifndef USEOPENMP
            for (Int_t i=1;i<numparts;i++)
            {
                temp =particle[i].Radius();
                if (temp<r[0])r[0]=temp;
                if (temp>r[2])r[2]=temp;
                r[1]+=temp*particle[i].GetMass();
                mtot+=particle[i].GetMass();
            }
#else
            Int_t i;
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:averad,mtot)
            for (i=1;i<numparts;i++)
            {
                averad+=particle[i].Radius()*particle[i].GetMass();
                mtot+=particle[i].GetMass();
            }
    }
            r[1]=averad;
            int nthreads;
    #pragma omp parallel 
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
            Double_t *mina=new Double_t[nthreads];
            Double_t *maxa=new Double_t[nthreads];
            for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=r[0];
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic,1) nowait
            for (i=1;i<numparts;i++)
            {
                temp =particle[i].Radius();
                if (temp<mina[omp_get_thread_num()])mina[omp_get_thread_num()]=temp;
                if (temp>maxa[omp_get_thread_num()])maxa[omp_get_thread_num()]=temp;
            }
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < r[0]) r[0] = mina[i];
            if (maxa[i] > r[2]) r[2] = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
#endif
        }
        r[1]/=mtot;
        return r;
    }


     //returns the mass weighted average distance from the CM
    Double_t System::AverageRadius(bool cmframe) const
    {
        Double_t r=0;
        Int_t i;
        if (!cmframe) {
            System S(*this);
            S.AdjustForCM();
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:r)
#endif
            for (i=0;i<numparts;i++) r+=S[i].GetMass()*S[i].Radius();
#ifdef USEOPENMP
    }
#endif
        }
        else {
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:r)
#endif
            for (i=0;i<numparts;i++) r+=particle[i].GetMass()*particle[i].Radius();
#ifdef USEOPENMP
    }
#endif
        }
        return r/TotalMass();
    }

    //returns the average velocity in the CM frame
    Double_t System::AverageSpeed(bool cmframe) const
    {
        Double_t v=0,mtot=0.;
        Int_t i;
        if (!cmframe) {
            System S(*this);
            S.AdjustForCMVel();
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:v,mtot)
#endif
            for (i=0;i<numparts;i++) {
                v+=S[i].GetMass()*S[i].AbsoluteVelocity();
                mtot+=S[i].GetMass();
            }
#ifdef USEOPENMP
    }
#endif
        }
        else {
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:v,mtot)
#endif
            for (i=0;i<numparts;i++){
                v+=particle[i].GetMass()*particle[i].AbsoluteVelocity();
                mtot+=particle[i].GetMass();
            }
#ifdef USEOPENMP
    }
#endif
        }
        return v/mtot;
    }

    //returns the average velocity in the CM frame
    Coordinate System::AverageVelocity(bool cmframe) const
    {
        Coordinate v(0,0,0);
        Double_t mtot=0.;
        Int_t i;
        if (!cmframe) {
            System S(*this);
            S.AdjustForCMVel();
#ifndef USEOPENMP
            for (i=0;i<numparts;i++) {
                for (int j=0;j<3;j++) v[j]+=S[i].GetMass()*S[i].GetVelocity(j);
                mtot+=S[i].GetMass();
            }
#else 
            Double_t vx=0.,vy=0.,vz=0.;
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:vx,vy,vz,mtot)
            for (i=0;i<numparts;i++) {
                vx+=S[i].GetMass()*S[i].Vx();
                vy+=S[i].GetMass()*S[i].Vy();
                vz+=S[i].GetMass()*S[i].Vz();
                mtot+=S[i].GetMass();
            }
    }
            v[0]=vx;v[1]=vy;v[2]=vz;
#endif
        }
        else {
#ifndef USEOPENMP
            for (i=0;i<numparts;i++){ 
                for (int j=0;j<3;j++) v[j]+=particle[i].GetMass()*particle[i].GetVelocity(j);
                mtot+=particle[i].GetMass();
            }
#else 
            Double_t vx=0.,vy=0.,vz=0.;
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:vx,vy,vz,mtot)
            for (i=0;i<numparts;i++) {
                vx+=particle[i].GetMass()*particle[i].Vx();
                vy+=particle[i].GetMass()*particle[i].Vy();
                vz+=particle[i].GetMass()*particle[i].Vz();
                mtot+=particle[i].GetMass();
            }
    }
            v[0]=vx;v[1]=vy;v[2]=vz;
#endif
        }
        return v*(1.0/mtot);
    }

    //returns the average radial velocity in the CM frame
    Double_t System::AverageRadialVelocity(bool cmframe) const
    {
        Double_t v=0,mtot=0.;
        Int_t i;
        if (!cmframe)
        {   
            System S(*this);
            S.AdjustForCMVel();
            S.AdjustForCM();
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:v,mtot)
#endif
            for (i=0;i<numparts;i++){
                v+=S[i].GetMass()*S[i].RadialVelocity();
                mtot+=S[i].GetMass();
            }
#ifdef USEOPENMP
    }
#endif
        }
        else {
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:v,mtot)
#endif
            for (i=0;i<numparts;i++){
                v+=particle[i].GetMass()*particle[i].RadialVelocity();
                mtot+=particle[i].GetMass();
            }
#ifdef USEOPENMP
    }
#endif
        }
        return v*(1.0/mtot);
    }

    //returns the max radial velocity in the CM frame
    Double_t System::MaxRadialVelocity(bool cmframe) const
    {
        Double_t v,temp;
        if (!cmframe)
        {
            System S(*this);
            S.AdjustForCMVel();
            S.AdjustForCM();
            v=S[0].RadialVelocity();
            for (Int_t i=1;i<numparts;i++)
            {
                temp=S[i].RadialVelocity();
                if(v<abs(temp))v=abs(temp);
            }
        }
        else {
            v=particle[0].RadialVelocity();
            for (Int_t i=1;i<numparts;i++)
            {
                temp=particle[i].RadialVelocity();
                if(v<abs(temp))v=abs(temp);
            }
        }
        return v;
    }

    //returns the max circular velocity in the CM frame
    Double_t System::MaxCircularVelocity(bool cmframe) const
    {
        Double_t v,temp;
        if (!cmframe)
        {
            System S(*this);
            S.AdjustForCMVel();
            S.AdjustForCM();
            v=S[0].CircularVelocity();
            for (Int_t i=1;i<numparts;i++)
            {
                temp=S[i].CircularVelocity();
                if(v<abs(temp))v=abs(temp);
            }
        }
        else {
            v=particle[0].CircularVelocity();
            for (Int_t i=1;i<numparts;i++)
            {
                temp=particle[i].CircularVelocity();
                if(v<abs(temp))v=abs(temp);
            }
        }
        return v;
    }

    //returns the average density
    Double_t System::AverageDensity() const
    {
        return TotalMass()/(4.0*3.14159/3.0*pow(MaxRadius(),3));
    }

    //returns the density in some radius R from some point O 
    Double_t System::RegionDensity(Double_t R, Double_t x, Double_t y, Double_t z) const
    {
        Double_t mass=0;
        for (Int_t i=0;i<numparts;i++)
        {
            if (sqrt(pow(particle[i].X()-x,2)+pow(particle[i].Y()-y,2)+pow(particle[i].Z()-z,2))<R)
                mass+=particle[i].GetMass();
        }
        return mass/(4.0*3.14159/3.0*pow(R,3));
    }

    //Sorting Functions
    void System::SortByID(){
        qsort(particle, numparts, sizeof(Particle), IDCompare);
    }
    void System::SortByRadius(){
        qsort(particle, numparts, sizeof(Particle), RadCompare);
    }
/*
    // Sort particles in system by radius
    void System::SortByRadius()
    {
        SortByRadius(0,numparts);
        //qsort(particle, numparts, sizeof(Particle), SystemCompare);
    }

    //is a quick sort. Eventually need to template this into 
    //a library of useful template sorts and other useful generic routines.
    void System::SortByRadius(int start, int n)
    {
        int pivot_index;
        int n1, n2;
        if (n > 1)
        {
            // partition inline
            Particle swap = particle[start + (n-1)/2]; // "middle" element
            particle[start + (n-1)/2] = particle[start];
            particle[start] = swap;
            Particle pivot = particle[start];
            Double_t pivot_rad = pivot.Radius();
            int too_big_index = start + 1;
            int too_small_index = start + n - 1;
            while (too_big_index <= too_small_index)
            {
                while (particle[too_big_index].Radius() <= pivot_rad && 
                       too_big_index < start + n)
                    too_big_index++; 
                while (particle[too_small_index].Radius() > pivot_rad)
                    too_small_index--;
                if (too_big_index < too_small_index)
                {
                    swap = particle[too_big_index];
                    particle[too_big_index] = particle[too_small_index];
                    particle[too_small_index] = swap;
                }
            }
            particle[start] = particle[too_small_index];
            particle[too_small_index] = pivot;
            pivot_index = too_small_index;
            // end inline partition
            n1 = pivot_index - start;
            n2 = n - n1 - 1;
            SortByRadius(start, n1);
            SortByRadius(pivot_index + 1, n2);
        }
    }
*/
    // Sort particles in system by distance to a postion
    void System::SortByDistance(Coordinate c)
    {
        this->SortByDistance(c,0,numparts);
        //qsort(particle, numparts, sizeof(Particle), SystemCompare);
    }

    //is a quick sort. Eventually need to template this into 
    //a library of useful template sorts and other useful generic routines.    
    void System::SortByDistance(Coordinate c, int start, int n)
    {
        int pivot_index;
        int n1, n2;

        if (n > 1)
        {
            // partition inline
            Particle swap = particle[start + (n-1)/2]; // "middle" element
            particle[start + (n-1)/2] = particle[start];
            particle[start] = swap;
            Particle pivot = particle[start];
            Double_t pivot_dist = sqrt(pow(pivot.X()-c[0],(Double_t)2.0)+pow(pivot.Y()-c[1],(Double_t)2.0)+pow(pivot.Z()-c[2],(Double_t)2.0));
            //cout <<c<<" "<<pivot_dist<<endl;
            int too_big_index = start + 1;
            int too_small_index = start + n - 1;

            while (too_big_index <= too_small_index)
            {
                Double_t dist;
                dist=sqrt(pow(particle[too_big_index].X()-c[0],(Double_t)2.0)+pow(particle[too_big_index].Y()-c[1],(Double_t)2.0)+pow(particle[too_big_index].Z()-c[2],(Double_t)2.0));
                while (dist <= pivot_dist && too_big_index < start + n)
                {
                    too_big_index++; 
                    dist=sqrt(pow(particle[too_big_index].X()-c[0],(Double_t)2.0)+pow(particle[too_big_index].Y()-c[1],(Double_t)2.0)+pow(particle[too_big_index].Z()-c[2],(Double_t)2.0));
                }

                dist=sqrt(pow(particle[too_small_index].X()-c[0],(Double_t)2.0)+pow(particle[too_small_index].Y()-c[1],(Double_t)2.0)+pow(particle[too_small_index].Z()-c[2],(Double_t)2.0));
                while (dist > pivot_dist)
                {
                    too_small_index--;
                    dist=sqrt(pow(particle[too_small_index].X()-c[0],(Double_t)2.0)+pow(particle[too_small_index].Y()-c[1],(Double_t)2.0)+pow(particle[too_small_index].Z()-c[2],(Double_t)2.0));
                }
                if (too_big_index < too_small_index)
                {
                    swap = particle[too_big_index];
                    particle[too_big_index] = particle[too_small_index];
                    particle[too_small_index] = swap;
                }
            }
    
            particle[start] = particle[too_small_index];
            particle[too_small_index] = pivot;
        
            pivot_index = too_small_index;
            // end inline partition
    
            n1 = pivot_index - start;
            n2 = n - n1 - 1;
        
            SortByDistance(c,start, n1);
            SortByDistance(c,pivot_index + 1, n2);
        }
    }

    void System::SortByEnergy(Double_t G, Double_t eps)
    {
        SortByEnergy(G,eps,0,numparts);
        //qsort(particle, numparts, sizeof(Particle), SystemCompare);
    }
    
    //is a quick sort. Eventually need to template this into 
    //a library of useful template sorts and other useful generic routines.    
    void System::SortByEnergy(Double_t G, Double_t eps, int start, int n)
    {
        int pivot_index;
        int n1, n2;
    
        if (n > 1)
        {
            // partition inline
            Particle swap = particle[start + (n-1)/2]; // "middle" element
            particle[start + (n-1)/2] = particle[start];
            particle[start] = swap;
            Particle pivot = particle[start];
            Double_t pivot_E = this->KineticEnergy(start)+G*(this->PotentialEnergy(start,eps));
            int too_big_index = start + 1;
            int too_small_index = start + n - 1;
            while (too_big_index <= too_small_index)
            {
                while ((this->KineticEnergy(too_big_index)+G*(this->PotentialEnergy(too_big_index,eps)) <= pivot_E) && 
                       (too_big_index < start + n))
                    too_big_index++; 
                while (this->KineticEnergy(too_small_index)+G*(this->PotentialEnergy(too_small_index,eps)) > pivot_E)
                    too_small_index--;
                if (too_big_index < too_small_index)
                {
                    swap = particle[too_big_index];
                    particle[too_big_index] = particle[too_small_index];
                    particle[too_small_index] = swap;
                }
            }
    
            particle[start] = particle[too_small_index];
            particle[too_small_index] = pivot;
        
            pivot_index = too_small_index;
            // end inline partition
    
            n1 = pivot_index - start;
            n2 = n - n1 - 1;
        
            SortByEnergy(G,eps, start, n1);
            SortByEnergy(G,eps, pivot_index + 1, n2);
        }
        //for (Int_t i=0;i<numparts;i++)
            //cout <<i<<" " << this->KineticEnergy(i)+G*(this->PotentialEnergy(i,eps)) ;
    }
    
    // Sort particles in system by radius
    void System::SortByDensity()
    {
        SortByDensity(0,numparts);
    }

    //is a quick sort. Eventually need to template this into 
    //a library of useful template sorts and other useful generic routines.    
    void System::SortByDensity(int start, int n)
    {
        int pivot_index;
        int n1, n2;
    
        if (n > 1)
        {
            // partition inline
            Particle swap = particle[start + (n-1)/2]; // "middle" element
            particle[start + (n-1)/2] = particle[start];
            particle[start] = swap;
            Particle pivot = particle[start];
            Double_t pivot_den = pivot.GetDensity();
            int too_big_index = start + 1;
            int too_small_index = start + n - 1;
    
            while (too_big_index <= too_small_index)
            {
                while (particle[too_big_index].GetDensity() <= pivot_den && 
                       too_big_index < start + n)
                    too_big_index++; 
                while (particle[too_small_index].GetDensity() > pivot_den)
                    too_small_index--;
                if (too_big_index < too_small_index)
                {
                    swap = particle[too_big_index];
                    particle[too_big_index] = particle[too_small_index];
                    particle[too_small_index] = swap;
                }
            }
    
            particle[start] = particle[too_small_index];
            particle[too_small_index] = pivot;
        
            pivot_index = too_small_index;
            // end inline partition
    
            n1 = pivot_index - start;
            n2 = n - n1 - 1;
        
            SortByDensity(start, n1);
            SortByDensity(pivot_index + 1, n2);
        }
    }

    //CM and reference frame functions
    
    // Subtract centre of mass coords ... needs testing, and
    // need to move to arrays rather than vectors
    Coord System::CM(Double_t tolerance) const
    {
        Coord Rcm = {0, 0, 0};
        Double_t TotalMass = 0.0;
        // Get centre of mass coords
        for (Int_t i = 0; i < numparts; i++)
        {
            for (int j = 0; j < 3; j++)
                Rcm.pos[j] += (particle[i].GetMass() * particle[i].GetPosition(j));
            TotalMass += particle[i].GetMass();
        }
        Double_t temp =1.0/TotalMass;
        for (int j = 0; j < 3; j++) Rcm.pos[j] *= temp;
        //if tolerance > 0 us interative technique.
        if (tolerance > 0.0)
        {
            // iterate to get a better c/m
            Coord Rcm_last = {Rcm.pos[0], Rcm.pos[1], Rcm.pos[2]};
            // find maximum radius of system centered on c/m
            Double_t x = particle[0].X() - Rcm.pos[0];
            Double_t y = particle[0].Y() - Rcm.pos[1];
            Double_t z = particle[0].Z() - Rcm.pos[2];
            Double_t rmax = sqrt(x * x + y * y + z * z);
            for (Int_t i = 1; i < numparts; i++)
            {
                x = particle[i].X() - Rcm.pos[0];
                y = particle[i].Y() - Rcm.pos[1];
                z = particle[i].Z() - Rcm.pos[2];
                Double_t test = sqrt(x * x + y * y + z * z);
                if (test > rmax) rmax = test;
            }
            Double_t ri = rmax;
            bool done = false;
            while (!done)
            {
                ri = 0.9 * ri;
                Double_t change = 1e32;
                while (change > tolerance)
                {
                    // find c/m of all particles within ri
                    Int_t Ninside = 0;
                    Double_t EncMass = 0.0;
                    for (int k = 0; k < 3; k++) Rcm.pos[k]=0.;
                    for (Int_t j = 0; j < numparts; j++)
                    {
                        x = particle[j].X() - Rcm_last.pos[0];
                        y = particle[j].Y() - Rcm_last.pos[1];
                        z = particle[j].Z() - Rcm_last.pos[2];
                        if (sqrt(x*x + y*y + z*z) <= ri)
                        {
                            for (int k = 0; k < 3; k++)
                                Rcm.pos[k] += particle[j].GetMass() *
                                          particle[j].GetPosition(k);
                            EncMass += particle[j].GetMass();
                            Ninside++;
                        }
                    }
                    for (int k = 0; k < 3; k++) Rcm.pos[k] /= EncMass;
                    // keep making radius smaller until there's
                    // less than 1% of the particle inside
                    if (EncMass < 0.01 * TotalMass) {done = true;break;}
                    // check for convergence
                    change = fabs(Rcm.pos[0] - Rcm_last.pos[0]);
                    for (int k = 1; k < 3; k++) if (fabs(Rcm.pos[k]-Rcm_last.pos[k])>change) change=fabs(Rcm.pos[k]-Rcm_last.pos[k]);
                    for (int k = 0; k < 3; k++) Rcm_last.pos[k] = Rcm.pos[k];
                }
            }
        }
        return Rcm;
    }

    // Subtract centre of mass coords 
    // need to optimize this
    void System::AdjustForCM(Double_t tolerance)
    {
        Coord Rcm=CM(tolerance);
        // Adjust particle positions
        for (Int_t i = 0; i < numparts; i++)
            particle[i].SetPosition(particle[i].GetPosition(0) - Rcm.pos[0],
                                    particle[i].GetPosition(1) - Rcm.pos[1],
                                    particle[i].GetPosition(2) - Rcm.pos[2]);
    }

    //returns velocity of centre of mass.
    Coord System::CMVel(Double_t tolerance) const
    {
        Coord Vcm = {0, 0, 0};
        Double_t TotalMass = 0.0;

        // Get centre of mass coords
        for (Int_t i = 0; i < numparts; i++)
        {
            for (int j = 0; j < 3; j++)
                Vcm.pos[j] += (particle[i].GetMass() * particle[i].GetVelocity(j));
            TotalMass += particle[i].GetMass();
        }
        Double_t temp =1.0/TotalMass;
        for (int j = 0; j < 3; j++)
            Vcm.pos[j]*= temp;
        return Vcm;
    }

     // Subtract centre of mass velocity... needs testing, and
    // need to move to arrays rather than vectors
    void System::AdjustForCMVel(Double_t tolerance)
    {
        Coord Vcm = CMVel();
        // Adjust particle positions
        for (Int_t i = 0; i < numparts; i++)
            particle[i].SetVelocity(particle[i].GetVelocity(0) - Vcm.pos[0],
                                    particle[i].GetVelocity(1) - Vcm.pos[1],
                                    particle[i].GetVelocity(2) - Vcm.pos[2]);
    }

    // Adjust particle positions
    // Subtract position (non-relativistic)
    void System::AdjustPosition(Coord R)
    {
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < numparts; i++)
            particle[i].SetPosition(particle[i].GetPosition(0) - R.pos[0],
                                    particle[i].GetPosition(1) - R.pos[1],
                                    particle[i].GetPosition(2) - R.pos[2]);
#ifdef USEOPENMP
    }
#endif

    }

    void System::AdjustPosition(Coordinate R)
    {
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < numparts; i++)
            particle[i].SetPosition(particle[i].GetPosition(0) - R[0],
                                    particle[i].GetPosition(1) - R[1],
                                    particle[i].GetPosition(2) - R[2]);
#ifdef USEOPENMP
    }
#endif
    }

    // Subtract velocity (non-relativistic)
    void System::AdjustVelocity(Coord V)
    {
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < numparts; i++)
            particle[i].SetVelocity(particle[i].GetVelocity(0) - V.pos[0],
                                    particle[i].GetVelocity(1) - V.pos[1],
                                    particle[i].GetVelocity(2) - V.pos[2]);
#ifdef USEOPENMP
    }
#endif
    }

    // Subtract velocity (non-relativistic)
    void System::AdjustVelocity(Coordinate V)
    {
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
        for (i = 0; i < numparts; i++)
            particle[i].SetVelocity(particle[i].GetVelocity(0) - V[0],
                                    particle[i].GetVelocity(1) - V[1],
                                    particle[i].GetVelocity(2) - V[2]);
#ifdef USEOPENMP
    }
#endif
    }

    //Angular Momentum
    Coord System::AngularMomentum(Coord o) const
    {
        Coord J = {0, 0, 0};
        Coord x = {0, 0, 0};
        Coord v = {0, 0, 0};
        for (Int_t i = 0; i < numparts; i++)
        {
            for (int j = 0; j < 3; j++){
                x.pos[j]=particle[i].GetPosition(j)-o.pos[0];
                v.pos[j]=particle[i].GetVelocity(j);
            }
            J.pos[0]+= particle[i].GetMass()*(x.pos[1]*v.pos[2]-x.pos[2]*v.pos[1]);
            J.pos[1]+= -particle[i].GetMass()*(x.pos[0]*v.pos[2]-x.pos[2]*v.pos[0]);
            J.pos[2]+= particle[i].GetMass()*(x.pos[0]*v.pos[1]-x.pos[1]*v.pos[0]);
        }
        return J;
    }
    Coord System::AngularMomentum(Coordinate o) const
    {
        Coord newo ={o[0],o[1],o[2]};
        Coord J = AngularMomentum(newo);
        return J;
    }
    //with no argument assumes center of mass
    Coord System::AngularMomentum() const
    {
        Coord o =this->CM();
        Coord J = AngularMomentum(o);
        return J;
    }

    //Energy Functions
    
    //Calculates the kinetic energy of the ith particle about the CM
    //Note that this KE has units of M*L^2/T^2
    Double_t System::KineticEnergy(Int_t i) const
    {
      //System S(*this);
      //S.AdjustForCMVel();

      Coord cmvel = CMVel();
      return 0.5*particle[i].GetMass()*(pow(particle[i].Vx()-cmvel.pos[0],2)
        +pow(particle[i].Vy()-cmvel.pos[1],2)+pow(particle[i].Vz()-cmvel.pos[2],2));
    }

    //Calculates the kinetic energy of the particles about the CM.
    //Note that this KE has units of M*L^2/T^2
    Double_t System::KineticEnergy() const
    {
      //note that making a copy of the system and adjusting its coordinates may use more memory but
      //might be more computationally efficient.
      //System S(*this);
      //S.AdjustForCMVel();

        Coord cmvel = CMVel();
        Double_t KE=0;
        for (Int_t i=0;i<numparts;i++)
            KE+=0.5*particle[i].GetMass()*(pow(particle[i].Vx()-cmvel.pos[0],2)
                +pow(particle[i].Vy()-cmvel.pos[1],2)+pow(particle[i].Vz()-cmvel.pos[2],2));
        //KE+=0.5*S[i].GetMass()*(pow(S[i].Vx(),2)+pow(S[i].Vy(),2)+pow(S[i].Vz(),2));
        return KE;

    }

    //Calculates the kinetic energy of the ith particle in a frame of reference.
    //Note that this KE has units of M*L^2/T^2
    Double_t System::KineticEnergy(Coord vel, Int_t i) const
    {
      return 0.5*particle[i].GetMass()*(pow(particle[i].Vx()-vel.pos[0],2)
        +pow(particle[i].Vy()-vel.pos[1],2)+pow(particle[i].Vz()-vel.pos[2],2));
    }

    //Calculates the kinetic energy of the particles in a frame of reference.
    //Note that this KE has units of M*L^2/T^2
    Double_t System::KineticEnergy(Coord vel) const
    {
        Double_t KE=0;
        for (Int_t i=0;i<numparts;i++)
        {
            KE+=0.5*particle[i].GetMass()*(pow(particle[i].Vx()-vel.pos[0],2)+
                pow(particle[i].Vy()-vel.pos[1],2)+pow(particle[i].Vz()-vel.pos[2],2));
        }
        return KE;

    }

    //Calculates the potential energy of the ith and jth particle
    //Note that to be allow any units, the potential energy is not multplied by G.
    //Thus V is in generic units of M^2/L.
    Double_t System::PotentialEnergy(Int_t i, Int_t j, Double_t eps) const
    {
        Double_t V=0;
        Double_t eps2=eps*eps;
            if (j!=i)
                V-=particle[i].GetMass()*particle[j].GetMass()*
                    1.0/pow((pow(particle[i].X()-particle[j].X(),(Double_t)2.0)
                    +pow(particle[i].Y()-particle[j].Y(),(Double_t)2.0)
                    +pow(particle[i].Z()-particle[j].Z(),(Double_t)2.0)+eps2),(Double_t)0.5);
            else
                V=0.0;
        return V;
    }
    
    //Calculates the potential energy of the ith particle in the gravtiational
    //Note that to be allow any units, the potential energy is not multplied by G.
    //Thus V is in generic units of M^2/L.
    Double_t System::PotentialEnergy(Int_t i, Double_t eps) const
    {
        Double_t V=0;
        Double_t eps2=eps*eps;
        for (Int_t j=0;j<numparts;j++)
        {
            //calculate potential between ith and jth particles. Note that G is not used
            if (j!=i)
                V-=particle[i].GetMass()*particle[j].GetMass()*
                    1.0/pow((pow(particle[i].X()-particle[j].X(),(Double_t)2.0)
                    +pow(particle[i].Y()-particle[j].Y(),(Double_t)2.0)
                    +pow(particle[i].Z()-particle[j].Z(),(Double_t)2.0)+eps2),(Double_t)0.5);
        }
        return V;
    }

    //Calculates the potential energy of the particles in the gravtiational
    //Note that to be allow any units, the potential energy is not multplied by G.
    //Thus V is in generic units of M^2/L.
    Double_t System::PotentialEnergy(Double_t eps, bool sphericalapprox) const
    {
        Double_t V=0;
        Double_t eps2=eps*eps;
        //the following code assumes some form of spherical symmetry.
        //its not that useful.
        if (sphericalapprox)
        {
            System S(*this);
            S.AdjustForCM();
            //S.AdjustForCMVel();
            S.SortByRadius();
            Double_t OutShell = 0.0;
            for (Int_t i = 0; i < numparts; i++)
            OutShell += S[i].GetMass() / sqrt(pow(S[i].Radius(),(Double_t)2.0) + eps2);

            Double_t EncMass = 0.0;
            for (Int_t i = 0; i < numparts; i++)
            {
                EncMass += S[i].GetMass();
                OutShell -= S[i].GetMass()/sqrt(pow(S[i].Radius(),(Double_t)2.0) + eps2);
                //V += (-1.0*EncMass/sqrt(pow(S[i].Radius(),2) + eps2) - OutShell);
                V += 0.5*S[i].GetMass()*(-1.0*EncMass/sqrt(pow(S[i].Radius(),(Double_t)2.0) + eps2) - OutShell);
            }
        }
        else
        {
            for (Int_t i=0;i<numparts;i++)
            {
                for (Int_t j=0;j<numparts;j++)
                {
                //calculate potential between ith and jth particles. the 0.5 is due to the Double_t
                //counting. Note that G is not used
                if (j!=i)
                    V-=0.5*particle[i].GetMass()*particle[j].GetMass()*
                        1.0/pow((pow(particle[i].X()-particle[j].X(),(Double_t)2.0)
                        +pow(particle[i].Y()-particle[j].Y(),(Double_t)2.0)+
                        pow(particle[i].Z()-particle[j].Z(),(Double_t)2.0)+eps2),(Double_t)0.5);
                }
            }
        }
        return V;
    }


    //Finds a particle and returns index. return -1 if not found.
    Int_t System::FindParticle(const Particle &part) const
    {
      Int_t i=0;
      while ((particle[i] != part)&&(i<numparts))
        i++;
      if (particle[i] != part) return -1;
      return i;
    }

    //Add or remove particles
    void System::AddParticle(Particle part)
    {
        Particle *testP=new Particle[numparts];
        for (Int_t i=0;i<numparts;i++)
            testP[i]=particle[i];
        delete[] particle;
        numparts++;
        particle = new Particle[numparts];
        for (Int_t i=0;i<numparts-1;i++)
            particle[i]=testP[i];
        particle[numparts-1]= part;
        delete []testP;
    }

    void System::RemoveParticle(Particle part)
    {
      if (numparts >0)
      {
        Int_t i=FindParticle(part);
        if (i==-1)
          cerr << "Can't find particle.\n";
        else
        {
            for (Int_t j=i;j<numparts-1;j++)
            {
                particle[j]=particle[j+1];
            }
          //delete &particle[numparts];
          numparts--;
          Particle *testP=new Particle[numparts];
          for (Int_t j=0;j<numparts;j++)
            testP[j]=particle[j];
          delete[] particle;
          particle = new Particle[numparts];
          for (Int_t j=0;j<numparts;j++)
            particle[j]=testP[j];
          delete[] testP;
        }
      }
      else
        cerr << "System is empty.\n";
    }

    void System::RemoveParticle(Int_t i)
    {
      if ((numparts >0) && (i>=0) && (i<numparts))
      {
          for (Int_t j=i;j<numparts-1;j++)
          {
            particle[j]=particle[j+1];
          }
          //delete &particle[numparts];
          numparts--;
          Particle *testP=new Particle[numparts];
          for (Int_t j=0;j<numparts;j++)
            testP[j]=particle[j];
          delete[] particle;
          particle = new Particle[numparts];
          for (Int_t j=0;j<numparts;j++)
            particle[j]=testP[j];
          delete[] testP;
      }
      else if (numparts==0)
        cerr << "System is empty.\n";
      else
        cerr << "Index out of range ("<<i<<").\n";
    }


    //Add system using vector between origins.
    void System::AddSystem(const System &s, const Double_t* o)
    {
        Particle *testP=new Particle[numparts];
        for (Int_t i=0;i<numparts;i++)
            testP[i]=particle[i];
        delete[] particle;
        for (Int_t i=0;i<numparts;i++)
            particle[i]=testP[i];

        particle = new Particle[numparts+s.numparts];

        for (Int_t i=numparts;i<(numparts+s.numparts);i++)
        {
          particle[i].SetMass(s.particle[i-numparts].GetMass());
          particle[i].SetPosition(s.particle[i-numparts].GetPosition(0)+o[0],
            s.particle[i-numparts].GetPosition(1)+o[1],
            s.particle[i-numparts].GetPosition(2)+o[2]);
          for (int j=0;j<3;j++)particle[i].SetVelocity(j,s.particle[i-numparts].GetVelocity(j));
        }
        numparts+=s.numparts;
        delete[] testP;
    }

    //Add system using vector between origins and relative velocity of origins.
    void System::AddSystem(const System &s, const Double_t* o, const Double_t* vo)
    {
        Particle *testP=new Particle[numparts];
        for (Int_t i=0;i<numparts;i++)
            testP[i]=particle[i];
        delete[] particle;
        for (Int_t i=0;i<numparts;i++)
            particle[i]=testP[i];

        particle = new Particle[numparts+s.numparts];

        for (Int_t i=numparts;i<(numparts+s.numparts);i++)
        {
          particle[i].SetMass(s.particle[i-numparts].GetMass());
          particle[i].SetPosition(s.particle[i-numparts].GetPosition(0)+o[0],
            s.particle[i-numparts].GetPosition(1)+o[1],
            s.particle[i-numparts].GetPosition(2)+o[2]);
          particle[i].SetVelocity(s.particle[i-numparts].GetVelocity(0)+vo[0],
            s.particle[i-numparts].GetVelocity(1)+vo[1],
            s.particle[i-numparts].GetVelocity(2)+vo[2]);
        }
        numparts+=s.numparts;
        delete[] testP;
    }

    // Remove all particles beyond a certain radius
    void System::ExtractSphere(Double_t radius)
    {
        //sort System by radius
        SortByRadius();
        Int_t i=0;
        while ((i<numparts)&&(particle[i].Radius()<radius)) i++;
        Particle *testP=new Particle[i+1];
        for (Int_t j=0;j<=i;j++)
            testP[j]=particle[j];
        numparts=i+1;
        delete[] particle;
        particle = new Particle[numparts];
        for (Int_t j=0;j<numparts;j++)
            particle[j]=testP[j];
        delete[] testP;
    }

    // Output system as binary data
    // ... to a file stream ...
    void System::Write(ofstream &F) const
    {
        //F.write((char *)&size, sizeof(int));
        F.write((char *)&numparts, sizeof(int));
        F.write((char *)&time, sizeof(Double_t));
        for (Int_t i=0;i<numparts;i++)
            particle[i].Write(F);
        //F.write((char *)particle, numparts * sizeof(Particle));
    }
    // ... to a file ...
    void System::Write(char *Filename) const
    {
        ofstream F(Filename, std::ios::out | std::ios::binary| std::ios::app);
        if (!F){
          cerr<<"Could not create file: "<<Filename<<endl;
          exit(8);
        }
        //F.write((char *)&size, sizeof(int));
        F.write((char *)&numparts, sizeof(int));
        F.write((char *)&time, sizeof(Double_t));
        for (Int_t i=0;i<numparts;i++)
            particle[i].Write(F);
        //F.write((char *)particle, numparts * sizeof(Particle));
        F.close();
    }
    // ... or to the screen.  This is c-style io.
    void System::Write(FILE *stream) const
    {
        //fwrite(&size, sizeof(int), 1, stream);
        fwrite(&numparts, sizeof(int), 1, stream);
        fwrite(&time, sizeof(Double_t), 1, stream);
        fwrite(particle, sizeof(Particle), numparts, stream);
    }


    //    OPERATORS
    System& System::operator=(const System &s)
    {
      delete [] particle;
      numparts=s.numparts;
      time = s.time;
      //size = s.size;
      particle = new Particle[numparts];
      for (Int_t i = 0; i<s.numparts; i++)
        particle[i]=s.particle[i];

      return *this;
    }

    //does size matter in equality?
    //equality does not care how the system is ordered.
    //this is etremely inefficient.
    bool System::operator==(const System &s) const
    {
      if ((numparts==s.numparts)&&(time==s.time))
      {
        Int_t j=0;
        do
        {
          if(FindParticle(s.particle[j])!=-1) j++;
          else return false;
        }while (j<s.numparts);
        return true;
      }
      return false;
    }

    //Have to think about how to warn if i is out of range.
    //What I've done here is not great.
    Particle &System::operator[](Int_t i)
    {
      if (i<numparts) return particle[i];
      else
      {
        cerr<<"Outside range, value must be < numparts ("<<numparts<<").";
        cerr<<"Returning last particle in array.";
        //return last member.
        return particle[numparts];
      }
    }

    Particle System::operator[](Int_t i) const
    {
      if (i<numparts) return particle[i];
      else
      {
        cerr<<"Outside range, value must be < numparts ("<<numparts<<").";
        cerr<<"Returning last particle in array.";
        //return last member.
        return particle[numparts];
      }
    }

    // Inserter: prints system as ascii data
    ostream &operator << (ostream &outs, const System &S)
    {
        outs << S.numparts << "\n";
        outs.setf(ios::scientific);
        outs.width(18); outs.precision(8);
        outs << S.time << "\n";
        outs.unsetf(ios::floatfield);

        if (S.numparts > 0)
        {
            for (Int_t i = 0; i < S.numparts; i++)
                outs << S.particle[i] << endl;
        }
        outs.unsetf(ios::floatfield);

        return outs;
    }

/*
    // Gas System
    //total mass in gas
    Double_t GasSystem::TotalGasMass() const
    {
      Double_t mt=0.;
      for (Int_t i=0;i<ngas;i++) mt+=gparticle[i].GetMass();
      return mt;
    }
    //returns mass weighted average temperature
    Double_t GasSystem::AverageTemp() const
    {
      Double_t avet=0, mt=0.;
      for (Int_t i=0;i<ngas;i++) {avet+=gparticle[i].GetTemp()*gparticle[i].GetMass();mt+=gparticle[i].GetMass();}
      return avet/mt;
    }
    Double_t GasSystem::MaxTemp() const
    {
      Double_t maxt=gparticle[0].GetTemp();
      for (Int_t i=1;i<ngas;i++) if (maxt<gparticle[i].GetTemp()) maxt=gparticle[i].GetTemp();
      return maxt;
    }
    */

    /* End of System */
}
