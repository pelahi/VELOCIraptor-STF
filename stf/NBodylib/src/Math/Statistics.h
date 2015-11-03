/*! \file Statistics.h
 *  \brief header file defining several statistical quantities
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream>
#include <complex>
#include <Precision.h>
#include <Coordinate.h>
#include <Matrix.h>
#include <GMatrix.h>
#include <cmath>

using namespace std;
namespace Math
{

    ///Average
    inline Double_t average(Int_t n, Double_t *x)
    {
        Double_t ave=0.;
        for (Int_t i=0;i<n;i++) ave+=x[i];
        return ave/(Double_t)n;
    }
    ///multidimensional average
    inline Double_t* average(Int_t parnum, Int_t datnum, Double_t **x)
    {
        Double_t *ave=new Double_t[parnum];
        for (Int_t i=0;i<parnum;i++) {
            ave[i]=0.;
            for (Int_t j=0;j<datnum;j++) ave[i]+=x[i][j];
            ave[i]/=(Double_t)datnum;
        }
        return ave;
    }
    ///Variance
    inline Double_t variance(Int_t n, Double_t xm, Double_t *x)
    {
        Double_t var=0.;
        for (Int_t i=0;i<n;i++) var+=(x[i]-xm)*(x[i]-xm);
        return var/(Double_t)(n-1.);
    }
    ///Variance without passing average
    inline Double_t variance(Int_t n, Double_t *x){ return variance(n, average(n,x), x);}
    ///Covariance matrix calculation
    inline GMatrix covariance(Int_t parnum, Int_t datnum, Double_t *xm,Double_t **x)
    {
        GMatrix covar(parnum,parnum);
        for (Int_t i=0;i<parnum;i++) 
            for (Int_t j=0;j<parnum;j++) {
                covar(i,j)=0.;
                for (Int_t k=0;k<datnum;k++) covar(i,j)+=(x[i][k]-xm[i])*(x[j][k]-xm[j]);
            covar(i,j)/=(Double_t)(datnum-1.0);
        }
        return covar;
    }
    ///Covariance matrix calculation without passing averages
    inline GMatrix covariance(Int_t parnum, Int_t datnum, Double_t **x)
    {
        return covariance(parnum,datnum,average(parnum,datnum,x),x);
    }
}

#endif
