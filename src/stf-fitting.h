/*! \file stf-fitting.h
 *  \brief this file contains inline function declarations for fitting the R distribution
 */

#ifndef STFFITTING_H
#define STFFITTING_H

///Math code
#include <NBodyMath.h>
using namespace Math;


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

inline int SkewGaussGSL(const gsl_vector *param, void *data, gsl_vector *f){
    gsl_fitting_data *fitdata = (gsl_fitting_data *)data;
    vector<double> params(fitdata->nallparams);
    for (int i=0;i<fitdata->nallparams;i++) params[i]=fitdata->allparams[i];
    for (int i=0;i<fitdata->npar;i++) params[fitdata->iparindex[i]] = gsl_vector_get(param, i);
    double x, y;
    for (int i=0; i<fitdata->n; i++) {
        x = fitdata ->x[i];
        if (x<=params[1]) y = params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]*params[3]));
        else y = params[0]*exp(-0.5*(x-params[1])*(x-params[1])/(params[2]));
        gsl_vector_set(f,i,y-fitdata->y[i]);
    }
    return GSL_SUCCESS;
}

inline int DiffSkewGaussAmpGSL(const gsl_vector *param, void *data, gsl_vector *f){
    gsl_fitting_data *fitdata = (gsl_fitting_data *)data;
    vector<double> params(fitdata->nallparams);
    for (int i=0;i<fitdata->nallparams;i++) params[i]=fitdata->allparams[i];
    for (int i=0;i<fitdata->npar;i++) params[fitdata->iparindex[i]] = gsl_vector_get(param, i);
    double ivar=1.0/params[2];
    double is2=1.0/params[3];
    double x, y, dx2;
    for (int i=0; i<fitdata->n; i++) {
        x = fitdata->x[i];
        dx2=(x-params[1])*(x-params[1]);
        if (x<=params[1]) y = exp(-0.5*dx2*ivar*is2);
        else y= exp(-0.5*dx2*ivar);
        gsl_vector_set(f,i,y);
    }
    return GSL_SUCCESS;
}

inline int DiffSkewGaussMeanGSL(const gsl_vector *param, void *data, gsl_vector *f){
    gsl_fitting_data *fitdata = (gsl_fitting_data *)data;
    vector<double> params(fitdata->nallparams);
    for (int i=0;i<fitdata->nallparams;i++) params[i]=fitdata->allparams[i];
    for (int i=0;i<fitdata->npar;i++) params[fitdata->iparindex[i]] = gsl_vector_get(param, i);
    double ivar=1.0/params[2];
    double is2=1.0/params[3];
    double x, y, dx2;
    for (int i=0; i<fitdata->n; i++) {
        x = fitdata->x[i];
        dx2=(x-params[1])*(x-params[1]);
        if (x<params[1]) y = params[0]*exp(-0.5*dx2*ivar*is2)*(x-params[1])*ivar*is2;
        else if (x==params[1]) y = 0.;
        else y = params[0]*exp(-0.5*dx2*ivar)*(x-params[1])*ivar;
        gsl_vector_set(f,i,y);
    }
    return GSL_SUCCESS;
}

inline int DiffSkewGaussVarGSL(const gsl_vector *param, void *data, gsl_vector *f){
    gsl_fitting_data *fitdata = (gsl_fitting_data *)data;
    vector<double> params(fitdata->nallparams);
    for (int i=0;i<fitdata->nallparams;i++) params[i]=fitdata->allparams[i];
    for (int i=0;i<fitdata->npar;i++) params[fitdata->iparindex[i]] = gsl_vector_get(param, i);
    double ivar=1.0/params[2];
    double is2=1.0/params[3];
    double x, y, dx2;
    for (int i=0; i<fitdata->n; i++) {
        x = fitdata->x[i];
        dx2=(x-params[1])*(x-params[1]);
        if (x<=params[1]) y = params[0]*exp(-0.5*dx2*ivar*is2)*dx2*ivar*ivar*is2*0.5;
        else y = params[0]*exp(-0.5*dx2*ivar)*dx2*ivar*ivar*0.5;
        gsl_vector_set(f,i,y);
    }
    return GSL_SUCCESS;
}

inline int DiffSkewGaussSkewGSL(const gsl_vector *param, void *data, gsl_vector *f){
    gsl_fitting_data *fitdata = (gsl_fitting_data *)data;
    vector<double> params(fitdata->nallparams);
    for (int i=0;i<fitdata->nallparams;i++) params[i]=fitdata->allparams[i];
    for (int i=0;i<fitdata->npar;i++) params[fitdata->iparindex[i]] = gsl_vector_get(param, i);
    double ivar=1.0/params[2];
    double is2=1.0/params[3];
    double x, y, dx2;
    for (int i=0; i<fitdata->n; i++) {
        x = fitdata->x[i];
        dx2=(x-params[1])*(x-params[1]);
        if (x<=params[1]) y = params[0]*exp(-0.5*dx2*ivar*is2)*dx2*ivar*is2*is2*0.5;
        else y = 0.0;
        gsl_vector_set(f,i,y);
    }
    return GSL_SUCCESS;
}

inline int DiffSkewGaussGSL(const gsl_vector *param, void *data, gsl_matrix *J){
    gsl_fitting_data *fitdata = (gsl_fitting_data *)data;
    gsl_vector * f = gsl_vector_alloc(fitdata->n);
    for (int iparam = 0; iparam < fitdata->npar; iparam++) {
        int activeparam = fitdata->iparindex[iparam];
        switch(activeparam)
        {
            case 0:
                DiffSkewGaussAmpGSL(param, data, f);
                break;
            case 1:
                DiffSkewGaussMeanGSL(param, data, f);
                break;
            case 2:
                DiffSkewGaussVarGSL(param, data, f);
                break;
            case 3:
                DiffSkewGaussSkewGSL(param, data, f);
                break;
        }
        for (int i=0;i<fitdata->n;i++) {
            gsl_matrix_set (J, i, iparam, gsl_vector_get(f,i));
        }
    }
    gsl_vector_free(f);
    return GSL_SUCCESS;
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
