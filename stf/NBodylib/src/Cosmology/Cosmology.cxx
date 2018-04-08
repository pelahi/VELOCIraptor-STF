#include <Cosmology.h>

using namespace std;
using namespace Math;

namespace Cosmology
{

    //scale factor function
    Double_t ScaleFactorFunc(const Double_t a, const Double_t om, const Double_t xla)
    {
        return 1/sqrt(1.+om*(1./a-1.)+xla*(pow(a,(Double_t)2.0)-1.));
    }

    Double_t HubbleFunc(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ola)
	{
		Double_t ok=1.0-om-ola;
		return Ho*sqrt(om/a/a/a+ola+ok/a/a);
	}
    Double_t HubbleFunc(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)
	{
		return Ho*sqrt(om/pow(a,Double_t(3.0*(1.0+alpha))+ok/a/a+ode/pow(a,(Double_t)(3.0*(1.0+wde)))));
	}

	//parameter at [10] specifices which HubbleFunction to call, 0 is standard, 1 is generalized.
	Double_t HubbleIntFunc(Double_t a, void* params)
	{
		Double_t * param=(Double_t*)params;
		if (param[10]<1)
		return 1.0/pow(a*HubbleFunc(a,1.0,param[1],param[2]),(Double_t)3.0);
		else return 1.0/pow(a*HubbleFunc(a,1.0,param[1],param[2],param[3],param[4],param[5]),(Double_t)3.0);
	}
	double HubbleIntFuncGSL(double a, void* params)
	{
//		double * param=(Double_t*)params;
		double * param=new double[20];
        for (int i=0;i<20;i++) param[i]=((Double_t*)params)[i];
		if (param[10]<1)
		return 1.0/pow(a*HubbleFunc(a,1.0,param[1],param[2]),3.0);
		else return 1.0/pow(a*HubbleFunc(a,1.0,param[1],param[2],param[3],param[4],param[5]),3.0);
	}
	double HubbleIntFuncGSLMonte(double a, size_t dim, void* params)
	{
//		double * param=(Double_t*)params;
		double * param=new double[20];
        for (int i=0;i<20;i++) param[i]=((Double_t*)params)[i];
		if (param[10]<1)
		return 1.0/pow(a*HubbleFunc(a,1.0,param[1],param[2]),3.0);
		else return 1.0/pow(a*HubbleFunc(a,1.0,param[1],param[2],param[3],param[4],param[5]),3.0);
	}

	Double_t aHIntFunc(Double_t a, void* params)
	{
		Double_t * param=(Double_t*)params;
		if (param[10]<1) return 1.0/(a*HubbleFunc(a,1.0,param[1],param[2]));
		else return 1.0/(a*HubbleFunc(a,1.0,param[1],param[2],param[3],param[4],param[5]));
	}
	double aHIntFuncGSL(double a, void* params)
	{
//		double * param=(Double_t*)params;
		double * param=new double[20];
        for (int i=0;i<20;i++) param[i]=((Double_t*)params)[i];
		if (param[10]<1) return 1.0/(a*HubbleFunc(a,1.0,param[1],param[2]));
		else return 1.0/(a*HubbleFunc(a,1.0,param[1],param[2],param[3],param[4],param[5]));
	}
	double aHIntFuncGSLMonte(double a, size_t dim, void* params)
	{
//		double * param=(Double_t*)params;
		double * param=new double[20];
        for (int i=0;i<20;i++) param[i]=((Double_t*)params)[i];
		if (param[10]<1) return 1.0/(a*HubbleFunc(a,1.0,param[1],param[2]));
		else return 1.0/(a*HubbleFunc(a,1.0,param[1],param[2],param[3],param[4],param[5]));
	}

    Double_t GrowthFunc(const Double_t a, const Double_t om, const Double_t ola)
	{
		Double_t HdivHo=HubbleFunc(a,1.0,om,ola);
		Double_t result;
		Double_t *param=new Double_t[20];
		param[0]=1.0;
		param[1]=om;
		param[2]=ola;
		param[10]=0;
		void *params=(void*)param;
		Double_t low[1],up[1];
		low[0]=0;
		up[0]=a;
		Double_t epsrel=1e-6;
		int numintervals=1000;

/*#ifdef GSL
	    gsl_monte_function gslFM1;
    	gslFM1.dim=1;
	    gslFM1.params=params;
#else*/
	    math_function F1;
    	F1.params=params;
//#endif

/*#ifdef GSL
    gslFM1.f=HubbleIntFuncGSL;
#else*/
    F1.function = HubbleIntFunc;
//#endif
	
//#ifdef ROMBERG
    	result = IntegrateRomberg(&F1,low[0],up[0],epsrel,numintervals,10);
//#elif CLOSED
//    	result = IntegrateClosed(&F1,low[0],up[0],numintervals);
//#elif GSL
//    	result = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
//#endif
		return 5.0/2.0*om*HdivHo*result;
	}

    Double_t GrowthFunc(const Double_t a, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)
	{
		Double_t HdivHo=HubbleFunc(a,1.0,om,ode,ok,alpha,wde);
		Double_t result;
		Double_t *param=new Double_t[20];
		param[0]=1.0;
		param[1]=om;
		param[2]=ode;
		param[3]=ok;
		param[4]=alpha;
		param[5]=wde;
		param[10]=1;
		void *params=(void*)param;
		Double_t low[1],up[1];
		low[0]=0;
		up[0]=a;
		Double_t epsrel=1e-6;
		int numintervals=1000;

/*#ifdef GSL
	    gsl_monte_function gslFM1;
    	gslFM1.dim=1;
	    gslFM1.params=params;
#else*/
	    math_function F1;
    	F1.params=params;
//#endif

/*#ifdef GSL
    gslFM1.f=HubbleIntFuncGSL;
#else*/
    F1.function = HubbleIntFunc;
//#endif
	
//#ifdef ROMBERG
    	result = IntegrateRomberg(&F1,low[0],up[0],epsrel,numintervals,10);
//#elif CLOSED
//    	result = IntegrateClosed(&F1,low[0],up[0],numintervals);
//#elif GSL
//    	result = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
//#endif
		return 5.0/2.0*om*HdivHo*result;
	}

    //Evaluates the time t (in units of t0) and Hubble parameter (units 1/t0)
    //given the expansion factor a (=1 at t=1)
    //for cosmological parameters omega0 and lambda0 (evaluated at t=1).
	Double_t Timet(Double_t a, const Double_t om, const Double_t ola)
	{
	    double h0t0,h0t;
	    Double_t *param=new Double_t[20];
		param[1]=om;
		param[2]=ola;
		param[10]=0;
		void*params= (void*) param;
		Double_t epsrel=1e-6;
		int numintervals=1000;
		params=(void*)param;

/*	    gsl_monte_function gslFM1;
  		gslFM1.dim=1;
    	gslFM1.params=params;
    	gsl_function gslF1;
	    gslF1.params=params;
*/	    //math_multidim_function FM1;
   		//FM1.ndim=1;
    	//FM1.params=params;
    	math_function F1;
	    F1.params=params;

    	/*if (logint)
		{
		    low[0]=log(low[0]);
    		up[0]=log(up[0]);
			gslFM1.f=aHfunclogintgslmonte;
	    	F1.function = aHfunclogint;
    	}*/
	    //else
    	{
//			gslFM1.f=aHIntFuncGSLMonte;
	    	F1.function = aHIntFunc;
    	}

		/*if (integrationtype==0){
	    	low[0]=0.0;up[0]=1.0;
    		h0t0 = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
	    	low[0]=0.0;up[0]=a;
    		h0t = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
		}*/
		//else if (integrationtype==1){
    		h0t0 = IntegrateRomberg(&F1,0.,1.0,epsrel,numintervals,10);
    		h0t = IntegrateRomberg(&F1,0.,a,epsrel,numintervals,10);
		//}
		/*else{
	    	h0t0 = IntegrateClosed(&F1,0.,1.0,numintervals);
	    	h0t = IntegrateClosed(&F1,0.,a,numintervals);
		}*/
    	return h0t/h0t0;
	}

    Double_t Timeh(const Double_t a, const Double_t om, const Double_t ola)
    {
	    double result;
	    Double_t *param=new Double_t[20];
		param[1]=om;
		param[2]=ola;
		param[10]=0;
		void*params= (void*) param;
		Double_t epsrel=1e-6;
		int numintervals=1000;
		params=(void*)param;

/*	    gsl_monte_function gslFM1;
  		gslFM1.dim=1;
    	gslFM1.params=params;
    	gsl_function gslF1;
	    gslF1.params=params;
*/    	math_function F1;
	    F1.params=params;

		//gslFM1.f=aHIntFuncGSLMonte;
    	F1.function = aHIntFunc;

/*		if (integrationtype==0){
	    	low[0]=0.0;up[0]=1.0;
    		result = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
		}*/
//		else if (integrationtype==1) 
		result = IntegrateRomberg(&F1,0.,1.0,epsrel,numintervals,10);
//		else result = IntegrateClosed(&F1,0.,1.0,numintervals);
        return result*sqrt(1.+om*(1./a-1.)+ola*(pow(a,2)-1.))/a;
    }
/*
    Double_t Timet(const Double_t a, const Double_t om, const Double_t ola)
    {
        Double_t h0t0=aintegral(0.,1.,om,ola);
        Double_t h0t=aintegral(0.,a,om,ola);
        return h0t/h0t0;
    }
    Double_t Timeh(const Double_t a, const Double_t om, const Double_t ola)
    {
        Double_t h0t0=aintegral(0.,1.,om,ola);
        Double_t h0t=aintegral(0.,a,om,ola);
        return h0t0*sqrt(1.+om*(1./a-1.)+ola*(pow(a,2)-1.))/a;
    }
*/
	Double_t WKR2(const Double_t k, const Double_t R){
    	Double_t kr=k*R;
    	Double_t kr6inv=1.0/(kr*kr*kr*kr*kr*kr);
    	Double_t s2=(sin(kr)-kr*cos(kr))*(sin(kr)-kr*cos(kr));
    	return 9.0*s2*kr6inv;
	}

	//power spectrums and their effective index
	//power-law
	Double_t PK(const Double_t k, const Double_t ns, const Double_t Amp) {return Amp*pow(k,ns);}

	//qBBKS
	Double_t qBBSK(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab){
		return k/(Omegam*h*h);
	}
	Double_t TBBKS(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab){
		Double_t q,gamma;
		//original BBKS q is 
		//q=k/(Omegam*h*h);
		//revised is 
		gamma=Omegam*h;
		if (gamma<=1e-4 &&ns>.1) gamma=Omegam*h*exp(-Omegab*(1.+1./Omegam));
        q=k/(gamma*h*h);
		return log(1+2.34*q)/(2.34*q)*pow(1+(3.89*q)+(16.1*q)*(16.1*q)+(5.4*q)*(5.4*q)*(5.4*q)+(6.71*q)*(6.71*q)*(6.71*q)*(6.71*q),-0.25);
	}
	//PBBSK here k is in Mpc^-1
	Double_t PBBKS(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Omegab){
		Double_t T=TBBKS(k,ns,Omegam,h,Omegab);
		return pow(k,ns)*Amp*T*T;
	}

	//----Based on Eisenstien and Hu 1998
	//q98
	Double_t qEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)
		{return k*Theta27*Theta27/Omegam/h/h;}
	//L98
	Double_t LEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)
		{return log(2.*exp(1.)+1.8*qEH98(k,Omegam,h,Theta27));}
	//C98
	Double_t CEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)
		{return 14.2+731./(1.+62.5*qEH98(k,Omegam,h,Theta27));}
	//T98
	Double_t TEH98(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)
		{return LEH98(k,Omegam,h,Theta27)/(LEH98(k,Omegam,h,Theta27)+CEH98(k,Omegam,h,Theta27)*pow(qEH98(k,Omegam,h,Theta27),(Double_t)2.0));}

	Double_t PEH98(const Double_t k, const Double_t ns, const Double_t deltaH, const Double_t clight, const Double_t Ho, const Double_t Omegam, const Double_t h, const Double_t Theta27){
    	Double_t temp=TEH98(k,Omegam,h,Theta27);
	    return 2.0*3.14159*3.14159*deltaH*deltaH*pow(clight/Ho,(Double_t)3.0+ns)*pow(k,ns)*temp*temp;
	}
	//P98 here k is in Mpc^-1, but just pass amplitude in
	Double_t PEH98(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27){
    	Double_t temp=TEH98(k,Omegam,h,Theta27);
    	return Amp*pow(k,ns)*temp*temp;
	}

	//P98 with WDM dampenning, here k is in Mpc^-1, but just pass amplitude in
	Double_t PWDMEH98(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd){
		Double_t damp=exp(-k*Rd/2.0-(k*Rd*k*Rd)/4.0);
    	return damp*damp*PEH98(k,ns,Amp,Omegam,h,Theta27);
	}

	Double_t LEH98diff(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)
		{return 1.8*Theta27*Theta27/Omegam/h/h/(2.0*exp(1.0)+1.8*qEH98(k,Omegam,h,Theta27));}
	//diff C98 by k
	Double_t CEH98diff(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)
		{return -45687.5/pow(1.+62.5*qEH98(k,Omegam,h,Theta27),2.)/Omegam/h/h*Theta27*Theta27;}
	//effective index in k units of Mpc^-1. neff =dlnP/dlnk
	Double_t neffEH98(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27)
	{
    	Double_t diff1=k*LEH98diff(k,Omegam,h,Theta27)/LEH98(k,Omegam,h,Theta27);
	    Double_t denom=1.0/(LEH98(k,Omegam,h,Theta27)+CEH98(k,Omegam,h,Theta27)*qEH98(k,Omegam,h,Theta27)*qEH98(k,Omegam,h,Theta27));
    	Double_t diff2=LEH98diff(k,Omegam,h,Theta27)+qEH98(k,Omegam,h,Theta27)*qEH98(k,Omegam,h,Theta27)*CEH98diff(k,Omegam,h,Theta27)
			+CEH98(k,Omegam,h,Theta27)*2.0*qEH98(k,Omegam,h,Theta27)*Theta27*Theta27/Omegam/h/h;
    	return ns+2.0*(diff1-diff2*denom*k);
	}

	Double_t neffWDMEH98(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)
	{
    	return neffEH98(k, ns, Omegam, h, Theta27)-(0.5*k*Rd)*(1.+k*Rd)*exp(-k*Rd/2.0-(k*Rd*k*Rd)/4.0);
	}

	//----Based on preprint of Green et al 04 for zeq>>z>zb and k>kb, zb=150 and kb~1e3 Mpc^-1.
	Double_t PGreen04(const Double_t k, const Double_t Amp, const Double_t wm){
		//assumes nu=2, c(nu)=3/2 b=-1.75 which is true only if f_b=Omega_b/Omega_m=0;
		Double_t keq=0.01*(wm/0.14);
		Double_t zeq=3371*(wm/0.14);
    	Double_t Dk=(1.0-2.0/3.0*(k*k/(1.7e6*1.7e6)))*exp(-(k*k/(1.7e6*1.7e6))-(k*k/(3.8e7*3.8e7)));
		Double_t lnterm=(log(k)-log(keq)-1.74)*(log(k)-log(keq)-1.74);
		Double_t zterm=(1.0+zeq)*(1.0+zeq);
		Double_t k3inv=1.0/(k*k*k);
		if (k>1e3)return 2.0*M_PI*k3inv*Amp*1.06*(1.5*1.5)*lnterm*Dk*zterm;
		else return 0.0;
	}

	//wrapper for integration of power spectrum *k^2;
	//simple power-law
	Double_t PKint(Double_t k, void *params)
	{
	    Double_t *param=(Double_t*) params;
    	return k*k*PK(k,param[0],param[1]);
	}
	double PKintgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
  	    return k*k*PK(k,param[0],param[1]);
	}
	Double_t PKlogint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PK(k,param[0],param[1]);
	}
	double PKlogintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
	    double *param=(double*) params;
    	return k*k*k*PK(k,param[0],param[1]);
	}

	//calls PBBKS
	Double_t PBBKSint(Double_t k, void *params)
	{
	    Double_t *param=(Double_t*) params;
    	return k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4]);
	}
	double PBBKSintgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
  	    return k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4]);
	}
	Double_t PBBKSlogint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4]);
	}
	double PBBKSlogintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
	    double *param=(double*) params;
    	return k*k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4]);
	}
    //calls amplitude PEH98
	Double_t PEH98int(Double_t k, void *params)
	{
	    Double_t *param=(Double_t*) params;
    	return k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4]);
	}
	double PEH98intgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
  	    return k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4]);
	}
	Double_t PEH98logint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4]);
	}
	double PEH98logintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
	    double *param=(double*) params;
    	return k*k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4]);
	}
    //PWDMEH98 wrappers
	Double_t PWDMEH98int(Double_t k, void *params)
	{
	    Double_t *param=(Double_t*) params;
    	return k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5]);
	}
	double PWDMEH98intgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
  	    return k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5]);
	}
	Double_t PWDMEH98logint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5]);
	}
	double PWDMEH98logintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
	    double *param=(double*) params;
    	return k*k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5]);
	}
	//Green power spectrum
	Double_t PGreen04int(Double_t k, void *params)
	{
    	Double_t *param=(Double_t*) params;
    	return k*k*PGreen04(k,param[0],param[1]);
	}
	double PGreen04intgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
	    double *param=(double*) params;
    	return k*k*PGreen04(k,param[0],param[1]);
	}
	Double_t PGreen04logint(Double_t logk, void *params)
	{
    	//cals amplitude PGreen04
    	Double_t k=exp(logk);
    	Double_t *param=(Double_t*) params;
    	return k*k*k*PGreen04(k,param[0],param[1]);
	}
	double PGreen04logintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
	    double *param=(double*) params;
    	return k*k*k*PGreen04(k,param[0],param[1]);
	}

	//integration of power spectrum, P(k)k^2dk*4pi/(2pi)^3
	Double_t IntegralkkPower(int powerspectype, void *params, int integrationtype, bool logint)
	{
	    double result;
    	double low[1],up[1];
	    int numintervals;
    	double epsrel;
	    double *param=(double*) params;
    	low[0]=param[10];
	    up[0]=param[11];
    	numintervals=(int)param[12];
	    epsrel=param[13];

	    gsl_monte_function gslFM1;
  		gslFM1.dim=1;
    	gslFM1.params=params;
	    //math_multidim_function FM1;
   		//FM1.ndim=1;
    	//FM1.params=params;
    	math_function F1;
	    F1.params=params;

    	if (logint)
		{
		    low[0]=log(low[0]);
    		up[0]=log(up[0]);
			//if (integrationtype==0){
			if (powerspectype==0) gslFM1.f=PEH98logintgslmonte;
			else if (powerspectype==1) gslFM1.f=PGreen04logintgslmonte;
			else if (powerspectype==2) gslFM1.f=PWDMEH98logintgslmonte;
			else if (powerspectype==3) gslFM1.f=PBBKSlogintgslmonte;
			else if (powerspectype==4) gslFM1.f=PKlogintgslmonte;
			//}
			//else {
	    	if (powerspectype==0) F1.function = PEH98logint;
	    	else if (powerspectype==1) F1.function = PGreen04logint;
	    	else if (powerspectype==2) F1.function = PWDMEH98logint;
	    	else if (powerspectype==3) F1.function = PBBKSlogint;
	    	else if (powerspectype==4) F1.function = PKlogint;
    //FM1.function=PGreen04logintmonte;
			//}
    	}
	    else
    	{
			//if (integrationtype==0){
			if (powerspectype==0) gslFM1.f=PEH98intgslmonte;
			else if (powerspectype==1) gslFM1.f=PGreen04intgslmonte;
			else if (powerspectype==2) gslFM1.f=PWDMEH98intgslmonte;
			else if (powerspectype==3) gslFM1.f=PBBKSintgslmonte;
			else if (powerspectype==4) gslFM1.f=PKintgslmonte;
			//}
			//else {
		    if (powerspectype==0) F1.function = PEH98int;
    		else if (powerspectype==1) F1.function = PGreen04int;
		    else if (powerspectype==2) F1.function = PWDMEH98int;
		    else if (powerspectype==3) F1.function = PBBKSint;
		    else if (powerspectype==4) F1.function = PKint;
		    //FM1.function=PGreen04intmonte;
			//}
    	}

		if (integrationtype==0)
    		result = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
		else if (integrationtype==1)
    		result = IntegrateRomberg(&F1,low[0],up[0],epsrel,numintervals,10);
		else
	    	result = IntegrateClosed(&F1,low[0],up[0],numintervals);

    	return result/(2.0*M_PI*M_PI);
	}

	//wrapper for integration of power spectrum *k^2*W(kR);
	//power law
	Double_t sigmaPKint(Double_t k, void *params)
	{
    	Double_t *param=(Double_t*) params;
    	return k*k*PK(k,param[0],param[1])*WKR2(k,param[15]);
	}
	double sigmaPKintgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
    	return k*k*PK(k,param[0],param[1])*WKR2(k,param[15]);
	}
	Double_t sigmaPKlogint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PK(k,param[0],param[1])*WKR2(k,param[15]);
	}
	double sigmaPKlogintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
    	double *param=(double*) params;
    	return k*k*k*PK(k,param[0],param[1])*WKR2(k,param[15]);
	}
	//PBBKS
	Double_t sigmaPBBKSint(Double_t k, void *params)
	{
    	Double_t *param=(Double_t*) params;
    	return k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	double sigmaPBBKSintgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
    	return k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	Double_t sigmaPBBKSlogint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	double sigmaPBBKSlogintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
    	double *param=(double*) params;
    	return k*k*k*PBBKS(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	//PEH98
	Double_t sigmaPEH98int(Double_t k, void *params)
	{
    	Double_t *param=(Double_t*) params;
    	return k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	double sigmaPEH98intgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
    	return k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	Double_t sigmaPEH98logint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	double sigmaPEH98logintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
    	double *param=(double*) params;
    	return k*k*k*PEH98(k,param[0],param[1],param[2],param[3],param[4])*WKR2(k,param[15]);
	}
	//PWDMEH98
	Double_t sigmaPWDMEH98int(Double_t k, void *params)
	{
    	Double_t *param=(Double_t*) params;
    	return k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5])*WKR2(k,param[15]);
	}
	double sigmaPWDMEH98intgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
    	return k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5])*WKR2(k,param[15]);
	}
	Double_t sigmaPWDMEH98logint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5])*WKR2(k,param[15]);
	}
	double sigmaPWDMEH98logintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
    	double *param=(double*) params;
    	return k*k*k*PWDMEH98(k,param[0],param[1],param[2],param[3],param[4],param[5])*WKR2(k,param[15]);
	}
	//Green 
	Double_t sigmaPGreen04int(Double_t k, void *params)
	{
    	Double_t *param=(Double_t*) params;
    	return k*k*PGreen04(k,param[0],param[1])*WKR2(k,param[15]);
	}
	double sigmaPGreen04intgslmonte(double *ka, size_t dim, void *params)
	{
    	double k=ka[0];
    	double *param=(double*) params;
    	return k*k*PGreen04(k,param[0],param[1])*WKR2(k,param[15]);
	}
	Double_t sigmaPGreen04logint(Double_t logk, void *params)
	{
    	Double_t k=exp(logk);
	    Double_t *param=(Double_t*) params;
    	return k*k*k*PGreen04(k,param[0],param[1])*WKR2(k,param[15]);
	}
	double sigmaPGreen04logintgslmonte(double *logk, size_t dim, void *params)
	{
    	double k=exp(logk[0]);
    	double *param=(double*) params;
    	return k*k*k*PGreen04(k,param[0],param[1])*WKR2(k,param[15]);
	}

	//integration of power spectrum, W^2(kR)P(k)k^2dk*4pi/(2pi)^3
	Double_t IntegralSigma(int powerspectype, Double_t R, void *params, int integrationtype, bool logint)
	{
	    double result;
    	double low[1],up[1];
	    int numintervals;
    	double epsrel;
	    double *param=(double*) params;
    	low[0]=param[10];
	    up[0]=param[11];
    	numintervals=(int)param[12];
	    epsrel=param[13];
		param[15]=R;
		params=(void*)param;

	    gsl_monte_function gslFM1;
  		gslFM1.dim=1;
    	gslFM1.params=params;
	    //math_multidim_function FM1;
   		//FM1.ndim=1;
    	//FM1.params=params;
    	math_function F1;
	    F1.params=params;

    	if (logint)
		{
		    low[0]=log(low[0]);
    		up[0]=log(up[0]);
			//if (integrationtype==0){
			if (powerspectype==0) gslFM1.f=sigmaPEH98logintgslmonte;
			else if (powerspectype==1) gslFM1.f=sigmaPGreen04logintgslmonte;
			else if (powerspectype==2) gslFM1.f=sigmaPWDMEH98logintgslmonte;
			else if (powerspectype==3) gslFM1.f=sigmaPBBKSlogintgslmonte;
			else if (powerspectype==4) gslFM1.f=sigmaPKlogintgslmonte;
			//}
			//else {
	    	if (powerspectype==0) F1.function = sigmaPEH98logint;
	    	else if (powerspectype==1) F1.function = sigmaPGreen04logint;
	    	else if (powerspectype==2) F1.function = sigmaPWDMEH98logint;
	    	else if (powerspectype==3) F1.function = sigmaPBBKSlogint;
	    	else if (powerspectype==4) F1.function = sigmaPKlogint;
    //FM1.function=PGreen04logintmonte;
			//}
    	}
	    else
    	{
			//if (integrationtype==0){
			if (powerspectype==0) gslFM1.f=sigmaPEH98intgslmonte;
			else if (powerspectype==1) gslFM1.f=sigmaPGreen04intgslmonte;
			else if (powerspectype==2) gslFM1.f=sigmaPWDMEH98intgslmonte;
			else if (powerspectype==3) gslFM1.f=sigmaPBBKSintgslmonte;
			else if (powerspectype==4) gslFM1.f=sigmaPKintgslmonte;
			//}
			//else {
		    if (powerspectype==0) F1.function = sigmaPEH98int;
    		else if (powerspectype==1) F1.function = sigmaPGreen04int;
		    else if (powerspectype==2) F1.function = sigmaPWDMEH98int;
		    else if (powerspectype==3) F1.function = sigmaPBBKSint;
		    else if (powerspectype==4) F1.function = sigmaPKint;
		    //FM1.function=PGreen04intmonte;
			//}
    	}

		if (integrationtype==0)
    		result = IntegrateVegasMonte(&gslFM1,low,up,numintervals);
		else if (integrationtype==1)
    		result = IntegrateRomberg(&F1,low[0],up[0],epsrel,numintervals,10);
		else
	    	result = IntegrateClosed(&F1,low[0],up[0],numintervals);

    	return result/(2.0*M_PI*M_PI);
	}

    //-------------------------------------------------------------------
    //  Does open-ended romberg integration of the scale factor integral from a1 to a2
    //  Based on a Numerical Recipes routine
    Double_t aintegral(const Double_t a1,const Double_t a2, const Double_t omega, const Double_t xlambda)
    {
        const Double_t eps=1.e-6;
        const int jmax=14;
        const int jmaxp=jmax+1; 
        const int k=5;
        Double_t *h=new Double_t[jmaxp+1];
        Double_t *s=new Double_t[jmaxp];
        Double_t x,del,ddel,sum,ss,dss;
        Double_t arange=a2-a1;
        int it;
        h[1]=1.;
        for (int j=1;j<=jmax;j++)
        {
            //midpoint evaluation
    	    if (j==1)
            {
               Double_t x=0.5*(a1+a2);
               s[j]=arange*1.0/(ScaleFactorFunc(x,omega,xlambda)*HubbleFunc(x,1.0,omega,xlambda));
               it=1;
            }
            else
            {
                int tnm=it;
                del=arange/(3.*tnm);
                ddel=2.*del;
                x=a1+0.5*del;
                sum=0.;
                for (int jt=0;jt<it;jt++)
                {
                    sum+=1.0/(ScaleFactorFunc(x,omega,xlambda)*HubbleFunc(x,1.0,omega,xlambda));//ScaleFactorFunc(x,omega,xlambda);
                    x+=ddel;
                    sum+=1.0/(ScaleFactorFunc(x,omega,xlambda)*HubbleFunc(x,1.0,omega,xlambda));//ScaleFactorFunc(x,omega,xlambda);
                    x+=del;
                }
                s[j]=(s[j-1]+arange*sum/tnm)/3.;
                it=3*it;
            }
            //extrapolate to zero step size
            if (j==k)
            {
                PolyInt(&h[j-k],&s[j-k],k,0.,ss,dss);
                if (fabs(dss)<=eps*fabs(ss)) return ss;
            }
            s[j+1]=s[j];
            h[j+1]=h[j]/9.;
        }
        cerr<<"aintegral:too many steps\n";
        exit(8);
    }
}
