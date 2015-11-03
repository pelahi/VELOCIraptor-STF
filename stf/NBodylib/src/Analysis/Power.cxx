
#include <iostream>
//#include <fftw3.h>
#include <NBody.h>
//#include <NBody/Analysis/Power.h>

namespace NBody
{
    void DensityGridNGP(System &S, Double_t *rho, int Ng, Double_t &L)
    {
    /*    int Ntot = Ng * Ng * Ng;
        for (int i = 0; i < Ntot; i++)
            rho[i] = 0.0;
        
        L = 0.0;
        for (int i = 0; i < S.NumParts(); i++)
        {
            if (S.Part(i).X() > L) L = S.Part(i).X();
            if (S.Part(i).Y() > L) L = S.Part(i).Y();
            if (S.Part(i).Z() > L) L = S.Part(i).Z();
        }
        std::cerr << "Box length is " << L << std::endl;
        // make box a little bigger so all particles are inside it
        L += 1e-10 * L;
        
        Double_t deltaX = L / Ng;
        for (int i = 0; i < S.NumParts(); i++)
        {
            int pos_x = int(S[i].X() / deltaX);
            int pos_y = int(S[i].Y() / deltaX);
            int pos_z = int(S[i].Z() / deltaX);
            
            rho[pos_z + Ng * (pos_y + Ng * pos_x)] +=
                    S[i].Mass() / (deltaX * deltaX * deltaX);
        }
    */}
    
    void CalcPowSpectrum(Double_t *rho, Double_t *p, int Ng, Double_t BoxLength)
    {
     /*   // FFT density to get delta_k's
        fftw_complex *dk = (fftw_complex *)
                           fftw_malloc(sizeof(fftw_complex) * Ng * Ng * 
                                                             (Ng/2 + 1));
        
        fftw_plan plan = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, rho, dk,
                                              FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        
        Double_t *delta = new Double_t[Ng * Ng * Ng];
        Double_t deltaK = 2.0 * M_PI / BoxLength;
        Double_t average = 0.0;
        int c = 0;
#define C (k + Ng*(j + Ng*i))
        for (int i = 0; i < Ng; i++)
        {
            for (int j = 0; j < Ng; j++)
            {
                for (int k = 0; k < Ng; k++)
                {
                    if (k < Ng/2 + 1)
                    {
                        delta[C] = (dk[c][0] * dk[c][0] + 
                                    dk[c][1] * dk[c][1]);
                        //std::cerr << "c = " << c << "\n";
                        c++;                        
                    }
                    else
                    {
                        //int cc = (Ng-k-1) + Ng*((Ng-j-1) + Ng*(Ng-i-1));
                        int cc = (Ng-k) + Ng*(j + Ng*i);
                        //std::cerr << "cc = " << cc << "\n";
                        delta[C] = delta[cc];
                    }                        
                    average += delta[C];
                }
            }
        }        
        Double_t volume = BoxLength * BoxLength * BoxLength;
        std::cerr << "Average <d_k^2> = " << average / volume << "\n";
        
        // average in spherical bins ...
        Double_t kx, ky, kz;
        for (int z = 0; z < Ng/2+1; z++)
        {
            Double_t wn = deltaK * z;
            
            p[z] = 0.0;
            c = 0;
            for (int i = 0; i < Ng; i++)
            {
                if (i < Ng/2+1) kx = deltaK * i;
                else kx = -deltaK * (Ng - i);
                for (int j = 0; j < Ng; j++)
                {
                    if (j < Ng/2+1) ky = deltaK * j;
                    else ky = -deltaK * (Ng - j);
                    for (int k = 0; k < Ng/2 + 1; k++)
                    {
                        kz = deltaK * k;
                        //std::cerr << "kx = " << kx << " ky = " << ky << " kz = " << kz << "\n";
                        Double_t kk = sqrt(kx*kx + ky*ky + kz*kz);
                        if (kk >= (wn - 0.5 * deltaK) && 
                            kk < (wn + 0.5 * deltaK))
                        {
                            p[z] += delta[c];
                        }
                        c++;
                    }
                }
            }
            Double_t volume = 2.0 / 3.0 * M_PI * (pow(wn + 0.5 * deltaK, 3.0) -
                            pow(wn - 0.5 * deltaK, 3.0));
            p[z] /= (volume * deltaK * deltaK * deltaK);
        }
        
        delete delta;
        fftw_free(dk);
    */}
}  
