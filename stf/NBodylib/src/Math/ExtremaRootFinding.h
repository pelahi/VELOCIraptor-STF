/*! \file ExtremaRootFinding.h
 *  \brief header file for finding roots and extrema of a function
    \todo not implemented yet
 */

#ifndef MINMAXROOT_H
#define MINMAXROOT_H

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

using namespace std;
namespace Math
{
    //// Brent bracket search for extrema
    //Double_t Brent();
}

#endif
