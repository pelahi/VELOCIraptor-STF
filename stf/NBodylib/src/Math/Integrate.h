/*! \file Integrate.h
 *  \brief header file for integrating functions or data.
 */

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <Precision.h>
#include <Function.h>
#include <Interpolate.h>
#include <Random.h>

using namespace std;
namespace Math
{
    /// Integrate function f(x) from x[lower] to x[upper].  Uses some
    /// algorithm I found somewhere, sometime.
    Double_t IntegrateData(const Double_t x[], const Double_t f[], const int lower, 
                         const int upper);

    /// Integrate function f(x) from x[lower] to x[upper].  Uses trapeziod
    /// algorithm.
    Double_t IntegrateTrap(const Double_t x[], const Double_t f[], const int a,
                         const int b);
    ///uses adaptive trapezoidal rule based on numerical recipes
    Double_t IntegrateTrapezoidal(const math_function * f, const Double_t a,
                         const Double_t b, const int n);
    Double_t IntegrateQTrap(const math_function * f, const Double_t a,
                         const Double_t b, const Double_t epsrel, int IMAX=16);
    ///uses simple trapezoidal
    Double_t IntegrateSimpleTrapezoidal(const math_function * f, const Double_t a,
                         const Double_t b, const int n);

    ///extended closed integration
    Double_t IntegrateClosed(const math_function * f, const Double_t a,
                         const Double_t b, const int n);

    ///simpsons 3/8 rule
    Double_t IntegrateSimpson(const math_function * f, const Double_t a,
                         const Double_t b);


    ///Romberg integration based on numerical recipes
    Double_t IntegrateRomberg(const math_function * f, const Double_t a,
                         const Double_t b, const Double_t epsrel, const int n, const int K);

    ///integration based vegas monte carlo integrator. can also return error and chisq
    ///vegas monte is either gsl or NR based vegas integrator.
    Double_t IntegrateVegasMonte(math_multidim_function * f, Double_t *a, Double_t *b, const int n, Double_t *error=0,Double_t*chisq=0,int interations=6);
    double IntegrateVegasMonte(gsl_monte_function * f, double *a, double *b, const int n, double *error=0,double *chisq=0);

    ///Romberg integration based on numerical recipes but uses vegas monte carlo integrator.
    Double_t IntegrateRombergMonte(math_multidim_function * f, Double_t *a,
                         Double_t *b, const Double_t epsrel, const int sampling, const int n, const int K, int interations=6);

    ///subroutine used in my nr vegas integrator.
    void rebin(Double_t rc, int nd, Double_t r[], Double_t xin[], Double_t xi[]);
    ///numerical recipes vegas integrator. for simplicity use IntegrateVegasMonte as interface
    void vegas(math_multidim_function *fxn, Double_t regn[], int ndim, int init,
    unsigned long ncall, int itmx, int nprn, long int idum, Double_t *tgral, Double_t *sd,
    Double_t *chi2a);
}

#endif
