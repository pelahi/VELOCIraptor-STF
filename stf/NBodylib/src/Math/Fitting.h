/*! \file Fitting.h
 *  \brief header file for fitting functions to data.
 */

#ifndef FITTING_H
#define FITTING_H

#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <Precision.h>
#include <Function.h>
#include <Matrix.h>
#include <GMatrix.h>
#include <Interpolate.h>
#include <ExtremaRootFinding.h>

using namespace std;
namespace Math
{
    /*! \brief Nonlinear fitting using Levenberg–Marquardt algorithm.

        Least Squares fitting of function f(x) from x[lower] to x[upper] to data y[lower]-y[upper]. Returns chi2 and altered params (which must be passed with initial guess)
        one can also pass a weight matrix (which is inverse of covariance matrix of data points), error on the relative change in chi2, and also a pointer that indicates whether 
        the parameters are held fixed. finally, if data binned dof reduced by 1 and can specify maximum number of iterations
        Note the function does not explicitly take into account errors on the x and y points save in the use of the weight matrix which is the inverse of the covariance matrix
        Furthermore, the covariance matrix for the parameters is determined using the confidence level assuming normaly distributed errors. If that is not the case, 
        the returned covariance matrix will not be correct and the degree to which it is correct cannot be correct for. Instead a more rigours determination of the confidence 
        region used to determine the covariance matrix, such as MC methods using the data as a sample of the underlying distribution and fitting resulting data numerous times
        must be done. The resulting distribution of parameters can then be fit itself using a function if necessary.
    */
    Double_t FitNonLinLS(const math_function fitfunc, const math_function *difffuncs, const int nparams, Double_t *params, GMatrix &covar,
                         const int npoints, const Double_t x[], const Double_t y[], GMatrix *Weight=NULL, Double_t error=1e-3, Double_t cl=0.95, int *fixparam=NULL, int binned=1, int maxiterations=1000, int iestimateerror=0);


    /*! \brief Given n points of data, determine optimal bin size and return binned data and bin centers. 
        
        This is done by dividing the data range into N bins of width Delta. Count the number of events k_i that enter the i'th bin.
        Calculate the mean and variance of the number of events as kmean=1/N *Sum k_i and kvar=1/N*Sum (k_i -kmean)^2
        Then calculate C(Delta)= (2*kmean-kvar)/Delta^2, a cost function
        Repeat while changing Delta and find Delta' that minimizes C(Delta). Return Delta'
        Based on Shimazaki H. and Shinomoto S., A method for selecting the bin size of a time histogram Neural Computation (2007) Vol. 19(6), 1503-1527 \n
        Can also provide weights for the data points.
    */
    Double_t OptimalBins(const int npoints, Double_t *points, Double_t xmin, Double_t xmax, Double_t *weights=NULL);
}

#endif
