#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <cmath>
//nbody code
#include <NBody.h>
#include <NBodyMath.h>
//gsl code
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

using namespace std;
using namespace Math;

namespace Cosmology
{
    //scale factor function
    Double_t ScaleFactorFunc(const Double_t a, const Double_t om, const Double_t xla);
    //{return 1/sqrt(1.+om*(1./a-1.)+xla*(pow(a,2.0)-1.));};

    //Hubble expansion factor, omegak is diff omega m omegala
    Double_t HubbleFunc(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ola);
    //Generalized Hubble expansion factor with omega_r=0 and omega_k specified
    Double_t HubbleFunc(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha=0.0, Double_t wde=-1.0);

	//wrapper for integration	
/*	Double_t HubbleIntFunc(Double_t a, void* params);
	double HubbleIntFuncGSL(double a, void* params);
	double HubbleIntFuncGSLMonte(double a, size_t dim, void* params);

	//wrapper for integration of 1/aH for time calculation.
	Double_t aHIntFunc(Double_t a, void* params);
	double aHIntFuncGSL(double a, void* params);
	double* aHIntFuncGSLMonte(double a, size_t dim, void* params);
*/
    //Evaluates the time t (in units of t0) and Hubble parameter (units 1/t0)
    //given the expansion factor a (=1 at t=1)
    //for cosmological parameters omega0 and lambda0 (evaluated at t=1).
    Double_t Timet(const Double_t a, const Double_t om, const Double_t ola);
    Double_t Timeh(const Double_t a, const Double_t om, const Double_t ola);

    //Growth factor
    Double_t GrowthFunc(const Double_t a, const Double_t om, const Double_t ola);
    Double_t GrowthFunc(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha=0.0, Double_t wde=-1.0);

	//window function
	Double_t WKR2(const Double_t k, const Double_t R);

	//Power Spectra
	Double_t PK(const Double_t k, const Double_t ns, const Double_t Amp);
	//----Based on Bardeen et al 1986 BBKS
	//qBBKS
	Double_t qBBSK(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab);
	Double_t TBBKS(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab);
	//PBBSK here k is in Mpc^-1
	Double_t PBBKS(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Omegab=0.);

	//----Based on Eisenstien and Hu 1998
	//q98
	Double_t qEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//L98
	Double_t LEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//C98
	Double_t CEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//T98
	Double_t TEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//P98 here k is in Mpc^-1. 
	Double_t PEH98(const Double_t k, const Double_t ns, const Double_t deltaH, const Double_t clight, const Double_t Ho, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//P98 here k is in Mpc^-1, but just pass amplitude in 
	Double_t PEH98(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//----P98 with Warm Dark Matter dampenning like term, here k is in Mpc^-1. 
	Double_t PWDMEH98(const Double_t k, const Double_t ns, const Double_t AMP, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd);

	//effective index in k units of Mpc^-1. neff =dlnP/dlnk
	//diff L98 by k
	Double_t LEH98diff(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//diff C98 by k
	Double_t CEH98diff(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	Double_t neffEH98(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27);
	//----same as above but with Warm Dark Matter dampenning like term, here k is in Mpc^-1. 
	Double_t neffWDMEH98(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd);


	//----Based on preprint of Green et al 04 for zeq>>z>zb and k>kb, zb=150 and kb~1e3 Mpc^-1.
	Double_t PGreen04(const Double_t k, const Double_t Amp, const Double_t wm);

	//----Based on data from a file not certain how to construct it yet (could use void * to store k and pk arrays.
	Double_t PFromFile(const Double_t k, void* params);

	// For integrals, void params has form of params[0]=ns;params[1]=Amp;params[2]=Omegam;params[3]=h;
	// Must eventually alter so that Omega_Lambda is not assumed to be 1-Omegam

	//integrate P(k)k^2dk*4pi/(2pi)^3
	Double_t IntegralkkPower(int powerspectype, void *params, int itype=0, bool logint=true);

	//integrate W^2(kR)P(k)k^2dk*4pi/(2pi)^3
	Double_t IntegralSigma(int powerspectype, Double_t R, void *params, int itype=0, bool logint=true);

    //-------------------------------------------------------------------
    //  Does open-ended romberg integration of the scale factor integral from a1 to a2
    //  Based on a Numerical Recipes routine
    Double_t aintegral(const Double_t a1,const Double_t a2, const Double_t omega, const Double_t xlambda);  
  
}      

#endif
