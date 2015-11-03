// InitCond.cpp copyright 2004 Joseph D. MacMillan

#include <iostream>
#include <cmath>
#include <string>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <fftw3.h>
#include <NBody.h>
#include <NBodyMath.h>
#include <InitCond.h>

using namespace Math;
using namespace NBody;
using namespace std;

namespace NBody
{
    System* CubicGrid(int N, Double_t BoxLength, Double_t Mtotal, Double_t time)
    {
        System *s = new System();
        
        Double_t deltaX = BoxLength / N;
        Double_t mass = Mtotal / (N * N * N);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    Double_t x = deltaX * i;
                    Double_t y = deltaX * j;
                    Double_t z = deltaX * k;
                    //Particle temp(mass, x, y, z, 0, 0, 0);
                    //s->AddParticle(temp);
                    s->AddParticle(Particle(mass, x, y, z, 0, 0, 0));
                }
            }
        }
        return s;
    }
    
    System* UniformSphere(int N, Double_t Radius, Double_t Mtotal, Double_t time)
    {
        System *s = new System();
        
        Double_t deltaX = Radius / N;
        Double_t mass = Mtotal / (N * N * N);
        for (int i = -N/2; i < N/2; i++)
        {
            for (int j = -N/2; j < N/2; j++)
            {
                for (int k = -N/2; k < N/2; k++)
                {
                    Double_t x = deltaX * i;
                    Double_t y = deltaX * j;
                    Double_t z = deltaX * k;
                    Double_t r=sqrt(x*x+y*y+z*z);
                    //Particle temp(mass, x, y, z, 0, 0, 0);
                    //add particle if it is within the radius
                    //if (r<Radius) s->AddParticle(temp);
                    if (r<Radius) s->AddParticle(Particle(mass, x, y, z, 0, 0, 0));
                }
            }
        }
        return s;
    }

/*    void CalcGaussDisp(int N, Double_t BoxLength, Double_t *Sx, Double_t *Sy,
                       Double_t *Sz, Cosmology &cosmo, unsigned int Seed)
    {
        // set up random numbers
        gsl_rng *RandomSeed = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(RandomSeed, Seed);
    
        // set k-grid
        fftw_complex *fx = (fftw_complex *)
                           fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
        fftw_complex *fy = (fftw_complex *)
                           fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
        fftw_complex *fz = (fftw_complex *)
                           fftw_malloc(sizeof(fftw_complex) * N * N * (N/2+1));
    
        Double_t deltaK = 2.0 * M_PI / BoxLength;
    
        // Gaussian random numbers
        Double_t a; 
        Double_t b;
    
        Double_t kx, ky, kz;
        int c = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N/2 + 1; k++)
                {
                    
                    kx = deltaK * i;
                    ky = deltaK * j;
                    
                    kz = deltaK * k;
                                
                    Double_t kk = sqrt(kx*kx + ky*ky + kz*kz);
            
                    if (kk == 0)
                    {
                        a = 0.0;
                        b = 0.0;
                    }
                    else 
                    {
                        a = sqrt(cosmo.PowerSpectrum(kk)) / kk / kk * 
                            gsl_ran_gaussian(RandomSeed, 1.0);
                        b = sqrt(cosmo.PowerSpectrum(kk)) / kk / kk * 
                            gsl_ran_gaussian(RandomSeed, 1.0);
                    }
        
                    fx[c][0] = 0.5 * kx * b;
                    fx[c][1] = 0.5 * kx * a;
            
                    fy[c][0] = 0.5 * ky * b;
                    fy[c][1] = 0.5 * ky * a;
                
                    fz[c][0] = 0.5 * kz * b;
                    fz[c][1] = 0.5 * kz * a;
                            
                    c++;
                }
            }
        }
       
        fftw_plan plan = fftw_plan_dft_c2r_3d(N, N, N, fx, Sx, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    
        plan = fftw_plan_dft_c2r_3d(N, N, N, fy, Sy, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    
        plan = fftw_plan_dft_c2r_3d(N, N, N, fz, Sz, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    
        fftw_free(fx);
        fftw_free(fy);
        fftw_free(fz);
    }
   
    void ApplyZelApprox(System &s, const Double_t *Sx, const Double_t *Sy,
                        const Double_t *Sz, Double_t BoxLength, Cosmology &cosmo,
                        bool GadgetVpec)
    {
        Double_t D = cosmo.GrowthRate();
        Double_t Ddot = cosmo.GrowthRateDot();
        Double_t average = 0.0;
        for (int i = 0; i < s.NumParts(); i++)
        {
            Double_t x = s[i].X() - D * Sx[i];
            Double_t y = s[i].Y() - D * Sy[i];
            Double_t z = s[i].Z() - D * Sz[i];
            
            average += fabs(D * Sx[i]);
            
            while (x > BoxLength) x -= BoxLength;
            while (x < 0.0) x += BoxLength;
            while (y > BoxLength) y -= BoxLength;
            while (y < 0.0) y += BoxLength;
            while (z > BoxLength) z -= BoxLength;
            while (z < 0.0) z += BoxLength;
            
            Double_t vx = -Ddot * Sx[i];
            Double_t vy = -Ddot * Sy[i];
            Double_t vz = -Ddot * Sz[i];
            
            if (GadgetVpec)
            {
                vx *= sqrt(1.0 + cosmo.Z());
                vy *= sqrt(1.0 + cosmo.Z());
                vz *= sqrt(1.0 + cosmo.Z());
            }
            
            s[i].SetPosition(x, y, z);
            s[i].SetVelocity(vx, vy, vz);
            
        }
        std::cerr << "ApplyZelApprox: Average x displacement is " <<
                     average/Double_t(s.NumParts()) << "\n";
    }*/
}            
