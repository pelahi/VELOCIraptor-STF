/*! \file stf-fitting.h
 *  \brief this file contains inline function declarations for fitting the R distribution
 */

#ifndef STFFITTING_H
#define STFFITTING_H

///\name Functions for Skew Gaussian distribution
///Here param[0] is amplitude, param[1] is mean, param[2] is variance and param[3] is s^2, the skew parameter
//@{
inline Double_t SkewGauss(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    if (x<=params[1]) return params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]*params[3]));
    else return params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]));
}
inline Double_t DiffSkewGaussAmp(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    Double_t dx2=(x-params[1])*(x-params[1]);
    Double_t ivar=1.0/params[2];
    Double_t is2=1.0/params[3];
    if (x<=params[1]) return exp(-0.5*dx2*ivar*is2);
    else return exp(-0.5*dx2*ivar);
}
inline Double_t DiffSkewGaussMean(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    Double_t dx2=(x-params[1])*(x-params[1]);
    Double_t ivar=1.0/params[2];
    Double_t is2=1.0/params[3];
    if (x<params[1]) return params[0]*exp(-0.5*dx2*ivar*is2)*(x-params[1])*ivar*is2;
    else if (x==params[1]) return 0.;
    else return params[0]*exp(-0.5*dx2*ivar)*(x-params[1])*ivar;
}
inline Double_t DiffSkewGaussVar(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    Double_t dx2=(x-params[1])*(x-params[1]);
    Double_t ivar=1.0/params[2];
    Double_t is2=1.0/params[3];
    if (x<=params[1]) return params[0]*exp(-0.5*dx2*ivar*is2)*dx2*ivar*ivar*is2*0.5;
    else return params[0]*exp(-0.5*dx2*ivar)*dx2*ivar*ivar*0.5;
}
inline Double_t DiffSkewGaussSkew(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    Double_t dx2=(x-params[1])*(x-params[1]);
    Double_t ivar=1.0/params[2];
    Double_t is2=1.0/params[3];
    if (x<=params[1]) return params[0]*exp(-0.5*dx2*ivar*is2)*dx2*ivar*is2*is2*0.5;
    else return 0.0;
}
//@}

///\name Functions for Gaussian distribution
///Here param[0] is amplitude, param[1] is mean, param[2] is variance
//@{
inline Double_t Gauss(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    return params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]));
}
inline Double_t DiffGaussAmp(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    return exp(-0.5*(x-params[1])*(x-params[1])/(params[2]));
}
inline Double_t DiffGaussMean(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    return params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]))*(x-params[1])/(params[2]);
}
inline Double_t DiffGaussVar(Double_t x, void *param){
    Double_t *params=(Double_t *)param;
    return params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]))*(x-params[1])*(x-params[1])/(params[2]*params[2])*0.5;
}
//@}
#endif
