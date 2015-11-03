/*! \file SpecialFunctions.h
 *  \brief header file defining several special functions
    This is a work in progress as generally can always call gsl for most of the functions.
 */

#ifndef SPECIALFUNCTIONS_H
#define SPECIALFUNCTIONS_H

#include <iostream>
#include <complex>
#include <Precision.h>
#include <cmath>
#include <gsl/gsl_sf.h>


using namespace std;
namespace Math
{

    //factorial
    int Factorial(int k);
    //Gamma function not yet idependently implemented (using just gsl but if want arbitrary precision must 
	//implement at some point.
	//complete gamma function
	//Double_t GammaFunc(Double_t a);
	//Calculates the normalized (regularized) upper incomplete gamma function (upper integral).
	//Double_t GammaFunc(Double_t a, Double_t x);

	//Again same issue with Gamma Functions, not yet implemented for arbitrary precision
    //calculates spherical harmonics
    //Double_t SphericalHarmonic(Double_t theta, Double_t phi, int l, int m);
    //just real component
    //Double_t SphericalHarmonicR(Double_t theta, Double_t phi, int l, int m);
    //just imaginary component
    //Double_t SphericalHarmonicI(Double_t theta, Double_t phi, int l, int m);

}

#endif
