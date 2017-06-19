/*! \file Fitting.cxx
 *  \brief subroutines for fitting functions to data.
 */

#include <Fitting.h>

using namespace std;
namespace Math
{
    ///\todo must alter so doesn't hang that often but I have no idea what could be taking so long, what hangs 
    Double_t FitNonLinLS(const math_function fitfunc, const math_function *difffuncs, const int nparams, Double_t *params, GMatrix &covar, 
                         const int npoints, const Double_t x[], const Double_t y[], GMatrix *W, Double_t error, Double_t cl, int *fixparam, int binned, int maxiit, int iestimateerror)
    {
        int iit=0,iflag,npar=0;
        int parlist[nparams];//store list of all unfixed parameters
        //first check if some parameters are held fixed
        if (fixparam!=NULL) {
            for (int i=0;i<nparams;i++) if (fixparam[i]==0) parlist[npar++]=i;
        }
        else {
            npar=nparams;
            for (int i=0;i<nparams;i++) parlist[i]=i;
        }
        int dof=npoints-npar-(binned==1);
        Double_t chi2,oldchi2,deltachi2,lambda=0.001,lambdabnd;
        GMatrix Jacobian(npoints,npar),JT(npoints,npar);
        GMatrix dY(npoints,1), dP(npar,1); //difference between data points and prediction, and resulting difference in parameters
        GMatrix a(npar), aa(npar), b(npar,1), ainv(npar);//temporary matrices
        iit=0;
        lambda=0.001;
        //check weights, if null initialize to identity matrix
        if (W==NULL) W=new GMatrix(npoints);
        do {
            ++iit;
            chi2=0;
            for (int i=0;i<npoints;i++){
                dY(i,0)=(y[i]-fitfunc.function(x[i],(void*)params));
                for (int j=0;j<npar;j++) Jacobian(i,j)=difffuncs[parlist[j]].function(x[i],(void*)params);
            }
            //check to see if getting nonsense
            for (int i=0;i<npoints;i++) {
              if (std::isnan(dY(i,0))) {return -1;}
            }
            chi2=(dY.Transpose()*(*W)*dY)(0,0);
            JT=Jacobian.Transpose();
            a=JT*(*W)*Jacobian;
            if (a.Trace()==0) {return -1;}
            lambdabnd=1.0/(a.Trace());
            //increase importance of diagonal elements of a with Levenberg-Marquardt parameter
            //aa=a.Diag();
            a=(a+a.Diag())*lambda;
            b=JT*(*W)*dY;
            //change in parameters is based on (J^T*W*J+lambda*diag(J^T*W*J))*dP=J^T*W*dY or a*dP=b, thus dP=a^-1*b
            ainv=a.InversewithPivot();
            covar=ainv;
            dP=ainv*b;
            //change parameters and check absolute change in values
            for (int j=0;j<npar;j++) params[parlist[j]]+=dP(j,0);
            //check to see if improved results
            oldchi2=chi2;chi2=0;
            for (int i=0;i<npoints;i++) dY(i,0)=(y[i]-fitfunc.function(x[i],(void*)params));
            chi2=(dY.Transpose()*(*W)*dY)(0,0);
            if (std::isnan(chi2)==false) {
                //if change from previous chi^2 minimum to current is below error tolorence and lambda is not large then
                //this is effectively a minimum
                deltachi2=(oldchi2-chi2)/chi2;
                deltachi2=deltachi2*deltachi2;
                if (chi2>oldchi2) {
                    lambda*=10.0;
                    for (int j=0;j<npar;j++) params[parlist[j]]-=dP(j,0);
                    iflag=0;
                }
                else if (deltachi2<error*error) iflag=1;
                else {
                    lambda/=10.0;
                    iflag=0;
                }
            }
            else {return -1;}
//                lambda/=10.0;
//                iflag=0;
//            }
            //if dampening parameter too big, fit failed
            if (lambda/MAXVALUE>0.1) {return -1;}
        }while (iflag==0&&iit<=maxiit);
        //if a fit was found meeting the expect variation requirements and is a local minimum (code does not guarantee finding global minimum)
        //then do one last fit with lambda =0; There appears to be an issue with whether lambda must be zero for appropriate fit.
        //and determine covariance matrix for parameters.
        //An estimate of the covariance matrix is the inverse of (JT*W*J) which should give you 1/sigma^2df/dPi*df/dPj
        //more accurate estimate requires determining the hessian matrix and inverting. Must make more robust call when determining covariance matrix
        //the parameters
        //if iterations have not converged print warning
        if (iit>maxiit) return -1;
        oldchi2=chi2;
        if (iestimateerror==1) {
            for (int i=0;i<npoints;i++)
                for (int j=0;j<npar;j++) Jacobian(i,j)=difffuncs[parlist[j]].function(x[i],(void*)params);
            JT=Jacobian.Transpose();
            a=JT*(*W)*Jacobian;
            b=JT*(*W)*dY;
            ainv=a.InversewithPivot();
            dP=ainv*b;
            for (int j=0;j<npar;j++) params[parlist[j]]+=dP(j,0);
            for (int i=0;i<npoints;i++) dY(i,0)=(y[i]-fitfunc.function(x[i],(void*)params));
            chi2=(dY.Transpose()*(*W)*dY)(0,0);
            if (chi2>oldchi2) {
                for (int j=0;j<npar;j++) params[parlist[j]]-=dP(j,0);
                chi2=oldchi2;
            }
            else {
                for (int i=0;i<npoints;i++)
                    for (int j=0;j<npar;j++) Jacobian(i,j)=difffuncs[parlist[j]].function(x[i],(void*)params);
                JT=Jacobian.Transpose();
                a=JT*(*W)*Jacobian;
                covar=a.InversewithPivot();
            }
            ///NOTE COVARIANCE ESTIMATOR STILL NOT WORKING!!!
            //also determine deltachi2 required for desired cl level which depends on dof
            deltachi2=gsl_cdf_chisq_Pinv(cl,dof);
        }
        //now alter covariance matrix to account for cl level
        //this has yet to be implemented but simply a matter of using deltachi2=dP*(covar^-1)*dP. Knowing this, can alter the values in the covarianc matrix to account for this change.
        return chi2;
    }


    Double_t OptimalBins(const int npoints, Double_t *points, Double_t xmin, Double_t xmax, Double_t *weights){
        //for bracketing to find mimimum, use very large number of bins and very small number
        Int_t nbins=max((Int_t)npoints/100.0,9.0)+1;
        Int_t nbins1=max((Int_t)npoints/100.0,9.0)+1,nbins2=5;
        Double_t *ki, *ki1,*ki2;
        Double_t delta1,delta2,deltamin,deltamid;
        Double_t kmean1,kvar1,kmean2,kvar2;
        Double_t Cdelta1,Cdelta2,dCddelta1,dCdelta2,Cdeltamin;
        
        //xbins=new Double_t[nbins];
        ki1=new Double_t[nbins1];
        ki2=new Double_t[nbins2];
        delta1=(xmax-xmin)/(Double_t)nbins1;
        delta2=(xmax-xmin)/(Double_t)nbins2;
        for (Int_t i=0;i<nbins1;i++) ki1[i]=0;
        for (Int_t i=0;i<nbins2;i++) ki2[i]=0;
        for (Int_t i=0;i<npoints;i++) {
            Int_t index1=(Int_t)(points[i]-xmin)/delta1;
            Int_t index2=(Int_t)(points[i]-xmin)/delta2;
            if (weights!=NULL) {
                ki1[index1]+=weights[i];
                ki2[index2]+=weights[i];
            }
            else {
                ki1[index1]++;
                ki2[index2]++;
            }
        }
        kmean1=kvar1=kmean2=kvar2=0.;
        for (Int_t i=0;i<nbins1;i++) kmean1+=ki1[i];
        kmean1/=(Double_t)nbins1;
        for (Int_t i=0;i<nbins1;i++) kvar1+=(ki1[i]-kmean1)*(ki1[i]-kmean1);
        kvar1/=(Double_t)nbins1;
        for (Int_t i=0;i<nbins2;i++) kmean2+=ki2[i];
        kmean2/=(Double_t)nbins2;
        for (Int_t i=0;i<nbins2;i++) kvar2+=(ki2[i]-kmean1)*(ki2[i]-kmean1);
        kvar2/=(Double_t)nbins2;
        //determine cost functions and initial bracket. Then begin minimum search using Brent's method
        Cdelta1=(2.0*kmean1-kvar1)/delta1;
        Cdelta2=(2.0*kmean2-kvar2)/delta2;
        Cdeltamin=min(Cdelta1,Cdelta2);
        if (Cdeltamin==Cdelta1)deltamin=delta1;
        else deltamin=delta2;
        
        //increase and decrease nbins and determine gradient of c with respect to delta 
        /*
        //now minimize cdelta based on initial estimate. Must pick direction 
        do {
            nbins*=fac
            ki=new Double_t[nbins];
            delta=(xmax-xmin)/(Double_t)nbins;
            for (Int_t i=0;i<nbins;i++) ki[i]=0;
            for (Int_t i=0;i<npoints;i++) {
                Int_t index=(Int_t)(points[i]-xmin)/delta;
                if (weights!=NULL) ki[index]+=weights[i];
                else ki[index]++;
            }
            kmean=kvar=0.;
            for (Int_t i=0;i<nbins;i++) kmean+=ki[i];
            kmean/=(Double_t)nbins;
            for (Int_t i=0;i<nbins;i++) kvar+=(ki[i]-mean)*(ki[i]-mean);
            kvar/=(Double_t)nbins;
            Cdelta=(2.0*kmean-kvar)/Delta;
            
        }while();*/
        return deltamin;
    }

}

