/*! \file Interpolate.h
 *  \brief header file for interpolation
 */

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <iostream> 
#include <cmath>
#include <NBodyMath.h>
#include <Precision.h>

using namespace std;
namespace Math
{
    ///polynomial interpolation to determine y value at a given x value baesd on data stored in x and y arrays
    void PolyInt(Double_t *xa,Double_t *ya,int n, Double_t x, Double_t &y, Double_t &dy);
}

#endif

