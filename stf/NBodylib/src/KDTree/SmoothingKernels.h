/*! \file SmoothingKernels.h
 *  \brief header file for inline smoothing kernel functions

    Several commonly used smoothing kernels. \n 
    Smoothing functions are normalized to unity using the volume enclosed by the function.
    For determining the normalization, one uses the corresponding volume element in a \a n-sphere, 
    where \a n is the dimensions of the sphere. That is for kernel \a U 
    \f[
        \int n C_n r^{n-1} dr U(r,h)=1,\quad C_n=\pi^{n/2}/\Gamma(n/2+1), 
    \f]
    where \f$ \Gamma \f$ is the gamma function.
    \todo Eventually alter so that kernal and dimension are separate since dimension just affects
    the normalization.

 */

#ifndef SMOOTHINGKERNELS_H
#define SMOOTHINGKERNELS_H

namespace NBody
{

    ///smoothing function structure
    struct smoothfunc {
        Double_t (* function) (Double_t r, Double_t h);
    };


    // Smoothing kernels.

    // Standard SPH smoothing kernel (this is a B-spline) (does not have normalization included as depends on dimension
    inline Double_t WSPH(Double_t r, Double_t h)
    {
        Double_t rh = r / h;
        if (rh <= 1.0) return (1.0-0.75*(2.0-rh)*(rh*rh));
        else if (rh>1.0&& rh <=2.0)return 0.25*(2.0-rh)*(2.0-rh)*(2.0-rh);
        else return 0;
    }
    // top hat kernel with R=2h
    inline Double_t WTH(Double_t r, Double_t h){
        Double_t rh=r/h*0.5;
        return (rh<=1.0);
    }
    // Gaussian smoothing kernel where sigma is 0.5*h
    inline Double_t WGauss(Double_t r, Double_t h)
    {
        Double_t rh = r / h*2.0;
        return exp(-0.5*pow(rh,(Double_t)2.0));
    }
    // Epanechnikov kernel (3/4*(1-x^2)) x <1, (note 3/4 just included in normalization)
    // so here x=r/2h ie, kernel smoothing scale is 2h
    inline Double_t WEpan(Double_t r, Double_t h){
        Double_t rh=r/h*0.5;
        if (rh<=1.0) return (1.0-rh*rh);
        else return 0.;
    }

}
#endif 
