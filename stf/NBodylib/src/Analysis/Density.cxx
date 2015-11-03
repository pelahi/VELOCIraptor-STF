
#include <iostream>
#include <NBody.h>
#include <Density.h>
#include <KDTree.h>
using namespace Math;

namespace NBody
{
    /*void DensityBinNorm(const System &S, int NumBins, Double_t *RBin, 
                        Double_t *DBin)
    {
        Double_t rmax = S.Rmax() + 1e-20;
        Double_t rmin = S.Rmin();
        
        Double_t deltaR = rmax / Double_t(NumBins);
        
        std::cerr << "Rmax = " << rmax << "\n";
        std::cerr << "Rmin = " << rmin << "\n";

        for (int i = 0; i < NumBins; i++)
        {
            DBin[i] = 0.0;
        }
        
        Double_t Mtot = 0.0;
        for (int i = 0; i < S.NumParts(); i++)
        {
            int pos = int(S[i].Radius() / deltaR);
            DBin[pos] += S[i].Mass();
            Mtot += S[i].Mass();
        }
        std::cerr << "Total mass is " << Mtot << std::endl;
        
        for (int i = 0; i < NumBins; i++)
        {
            // RBin is the average radius of each shell
            RBin[i] = 0.5 * deltaR + deltaR * i;
            Double_t R = deltaR * i;
            Double_t volume = 4.0 / 3.0 * M_PI * 
                           (pow(R + deltaR, 3.0) - pow(R, 3.0));
            DBin[i] /= volume;
        }
    }
    
    void DensityBinLog(const System &S, int NumBins, Double_t *RBin, 
                       Double_t *DBin)
    {
        Double_t rmax = S.Rmax() + 1e-20;
        Double_t rmin = S.Rmin();
        
        Double_t deltaR = log10(rmax / rmin) / Double_t(NumBins);
        
        std::cerr << "Rmax = " << rmax << "\n";
        std::cerr << "Rmin = " << rmin << "\n";
        
        for (int i = 0; i < NumBins; i++)
        {
            DBin[i] = 0.0;
        }
        
        for (int i = 0; i < S.NumParts(); i++)
        {
            int pos = int((log10(S[i].Radius()) - log10(rmin)) / deltaR);
            DBin[pos] += S[i].Mass();
        }
        
        for (int i = 0; i < NumBins; i++)
        {
            RBin[i] = pow(10.0, log10(rmin) + 0.5 * deltaR + deltaR * i);
            Double_t R = pow(10.0, deltaR * i) * rmin;
            Double_t Rplus1 = pow(10.0, deltaR) * R;
            Double_t volume = 4.0 / 3.0 * M_PI * 
                           (pow(Rplus1, 3.0) - pow(R, 3.0));
            DBin[i] /= volume;
        }   
    }
    
    void DensityBinEqualMass(System &S, int NumBins, Double_t *RBin, 
                             Double_t *DBin)
    {
        // how many particles in each bin?
        int N = int(ceil(Double_t(S.NumParts()) / Double_t(NumBins)));
        std::cerr << "Putting " << N << " particles into each bin\n";
        
        S.SortByRadius();
        
        Double_t rmax = S[S.NumParts() - 1].Radius() + 1e-10;
        Double_t rmin = S[0].Radius();
        std::cerr << "Rmax = " << rmax << "\n";
        std::cerr << "Rmin = " << rmin << "\n";
        
        int pos = 0;
        Double_t m = 0.0;
        Double_t r;
        Double_t r_last = rmin;
        for (int i = 0; i < S.NumParts(); i++)
        {
            m += S[i].Mass();
            if (((i+1) % N) == 0)
            {
                r = S[i].Radius();
                RBin[pos] = r_last + 0.5 * (r - r_last);
                Double_t volume = 4.0 / 3.0 * M_PI * (pow(r, 3.0) - 
                                pow(r_last, 3.0));
                DBin[pos] = m / volume;
                m = 0.0;
                r_last = r;
                pos++;
            }
        }
        
        if (pos == NumBins - 1)
        {
            r = rmax;
            RBin[pos] = r_last + 0.5 * (r - r_last);
            Double_t volume =  4.0 / 3.0 * M_PI * (pow(r, 3.0) - 
                             pow(r_last, 3.0));
            DBin[pos] = m / volume;
        }
    }
    
    void DensityBinLogEllip(const System &S, int NumBins, Double_t *RBin, 
                            Double_t *DBin, Double_t q, Double_t s)
    {
        Double_t rmax = S.RmaxEllip(q, s) + 1e-20;
        Double_t rmin = S.RminEllip(q, s);
        
        Double_t deltaR = log10(rmax / rmin) / Double_t(NumBins);
        
        std::cerr << "Rmax = " << rmax << "\n";
        std::cerr << "Rmin = " << rmin << "\n";
        
        for (int i = 0; i < NumBins; i++)
        {
            DBin[i] = 0.0;
        }
        
        for (int i = 0; i < S.NumParts(); i++)
        {
            int pos = int((log10(S[i].RadiusEllip(q, s)) - 
                           log10(rmin)) / deltaR);
            DBin[pos] += S[i].Mass();
        }
        
        for (int i = 0; i < NumBins; i++)
        {
            RBin[i] = pow(10.0, log10(rmin) + 0.5 * deltaR + deltaR * i);
            Double_t R = pow(10.0, deltaR * i) * rmin;
            Double_t Rplus1 = pow(10.0, deltaR) * R;
            Double_t volume = 4.0 / 3.0 * M_PI * 
                           (pow(Rplus1, 3.0) - pow(R, 3.0));
            DBin[i] /= volume;
        }   
    }*/

    void DensitySmooth(System &S, int BucketSize, int NumNearest)
    {
        // Create a tree from S
        KDTree *tree = new KDTree(S, BucketSize);
        // Calculate density at each particle, storing it in the particle data.
        tree->CalcDensity(NumNearest);
        delete tree;
    }

	/*void DensitySmoothBall(System &S, Double_t ballsize)
	{
		Int_t numparts =S.GetNumParts();
        for (Int_t i = 0; i < numparts; i++) S[i].SetDensity(0);

        Double_t furthest = INFINITY;

		//for each particle 
        for (Int_t i = 0; i < numparts; i++)
        {

		}

	}*/
}  
