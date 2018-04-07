/*! \file Integrate.cxx
 *  \brief integration subroutines 
 */

#include <cmath>
#include <Integrate.h>

using namespace std;
namespace Math
{
    Double_t IntegrateData(const Double_t x[], const Double_t f[], const int lower, 
                         const int upper)
    {
        Double_t  value = 0.0, cof;
        Double_t  first1, first2, second1, second2, fac1, fac2;
        int     i;

        first1 = (f[lower + 1] - f[lower]) / (x[lower + 1] - x[lower]);

        // get rid of warnings, at least
        first2 = second1 = second2 = 0.0;

        for (i = lower + 1; i <= upper; i++)
        {
            if (i < upper)
            {
                first2 = (f[i + 1] - f[i]) / (x[i + 1] - x[i]);
                second2 = (first2 - first1) / (x[i + 1] - x[i - 1]);
            }

            if (i == upper || i == lower + 1)
                second1 = second2;

            cof = x[i] - x[i - 1];
            fac1 = .5 * (f[i - 1] + f[i]);
            fac2 = .5 * (second1 + second2);

            value += cof * (fac1 - cof * cof * fac2 / 6.0);

            first1 = first2;
            second1 = second2;
        }

        return (value);
    }

    Double_t IntegrateTrap(const Double_t x[], const Double_t f[], const int a,
                         const int b)
    {
        int     i;
        Double_t  Tn = 0.0;

        if (a == b)
        {
            if (a == 0)
                return f[a] * (x[a + 1] - x[a]);
            else
                return f[a] * (x[a] - x[a - 1]);
        }

        for (i = a; i < b; ++i)
        {
            Tn += 0.5 * (x[i + 1] - x[i]) * (f[i + 1] + f[i]);
        }

        return Tn;
    }

    //based on Numeric Recipes
    Double_t IntegrateTrapezoidal(const math_function * f, const Double_t a,
                         const Double_t b, const int n)
    {
        int it,j;
        Double_t  sum,s, del,x,tnm;
        //cerr<<n<<" "<<a<<" "<<((Double_t *)f->params)[0]<<" "<<f->function(a,f->params)<<endl;
        if (n<=1) return (0.5*(b-a)*(f->function(a,f->params)+f->function(b,f->params)));
        /*else
        {
            si=f->function(a,f->params)+f->function(b,f->params))
            tnm=(Double_t)it;
          	del=(b-a)/(Double_t)n;
            for (sum=si,j=1;j<=n;j++){
				x=a+0.5*del;
				sum+=f->function(x,f->params);
            s=0.5*(s+(b-a)*sum/tnm);
            return 2.0*s;
        }*/
        else
        {
            for (j=1,it=1;j<n-1;j++)it<<=1;
            tnm=(Double_t)it;
            del=(b-a)/tnm;
            x=a+0.5*del;
            for (sum=0.0,j=1;j<=it;j++,x+=del)sum+=f->function(x,f->params);
            s=0.5*(s+(b-a)*sum/tnm);
            return 2.0*s;
        }
    }

    //based on Numeric Recipes
    Double_t IntegrateQTrap(const math_function * f, const Double_t a,
                         const Double_t b, const Double_t epsrel, int IMAX)
    {
        int     i;
        Double_t  s,olds;

        olds=-1.0e30;
        //cerr<<IMAX<<endl;
        //any number that is unlikely to be average of end points
        for (i=1;i<IMAX;i++)
        {
            //s=IntegrateTrapezoidal(f,a,b,i);
            s=IntegrateSimpleTrapezoidal(f,a,b,(int)pow(2.0,i-1));
            //cerr<<s<<endl;
            //if (fabs(s-olds)<epsrel*fabs(olds)) return s;
            //if (fabsl(s-olds)<epsrel*fabsl(olds)) return s;
            if ((s-olds)*(s-olds)<epsrel*epsrel*(olds*olds)) return s;
            if (s==0.0 && olds ==0.0 && i>IMAX/2) return s;
            olds=s;
        }
        cerr<<"No convergence\n";
        return s;
    }


    Double_t IntegrateSimpleTrapezoidal(const math_function * f, const Double_t a,
                         const Double_t b, const int n)
    {
        Double_t  sum,del,x,arg;
        del=(b-a)/(Double_t)(n-1);
        sum=0.0;
        sum+=0.5*del*(f->function(a,f->params)+f->function(b,f->params));
        for (int j=1;j<n-1;j++)
        {
            x=a+del*j;
            arg=f->function(x,f->params);
          	sum+=del*arg;
        }
        return sum;
    }

    Double_t IntegrateClosed(const math_function * f, const Double_t a,
                         const Double_t b, const int n)
    {
        Double_t  sum,del,x,arg;
        //n must be => 6)
        if (n<6) del=(b-a)/(Double_t)(5);
        else del=(b-a)/(Double_t)(n-1);
        sum=0.0;
        //first 3 points
        sum+=del*(3.0/8.0*f->function(a,f->params)+7.0/6.0*f->function(a+del,f->params)+23.0/24.0*f->function(a+2.0*del,f->params));

        //last 3 points
        sum+=del*(3.0/8.0*f->function(b,f->params)+7.0/6.0*f->function(b-del,f->params)+23.0/24.0*f->function(b-2.0*del,f->params));
        //rest
        for (int j=3;j<n-3;j++)
        {
            x=a+del*j;
            arg=f->function(x,f->params);
          	sum+=del*arg;
        }
        return sum;
    }

    Double_t IntegrateSimpson(const math_function * f, const Double_t a,
                         const Double_t b)
    {
        Double_t  sum,del;
        del=(b-a)/(Double_t)(3);
        sum=0.0;
        sum=del*(3.0/8.0*f->function(a,f->params)+9.0/8.0*f->function(a+del,f->params)+9.0/8.0*f->function(a+2.0*del,f->params)+3.0/8.0*f->function(b,f->params));
        return sum;
    }

    //
    Double_t IntegrateRomberg(const math_function * f, const Double_t a,
                         const Double_t b, const Double_t epsrel, const int n, const int nextrap)
    {
        Double_t ss,dss;
        //stores trapezoidal approximations and stepsizes
        Double_t *s=new Double_t[n+1];
        Double_t *h=new Double_t[n+2];
        //note that 
        h[0]=1.0;
        for (int j=0;j<n;j++)
        {
            s[j]=IntegrateTrapezoidal(f,a,b,j+1);
            //s[j]=IntegrateSimpleTrapezoidal(f,a,b,(int)pow(2.0,j+1));
            //s[j]=IntegrateClosed(f,a,b,j+6);
            if (j>=nextrap-1){
                    PolyInt(&h[j-nextrap+1], &s[j-nextrap+1], nextrap, 0.0, ss, dss);
                    //if (fabs(dss)<=epsrel*fabs(ss))return ss;
                    if (dss*dss<=epsrel*epsrel*ss*ss)return ss;
                    //if ((dss*dss)<=epsrel*epsrel*(ss*ss))return ss;
            }
            h[j+1]=0.25*h[j];
            //factor of 0.25 even though stepsize only decreaseed by 0.5 makes polynomial in extrapolate decrease in step size 
        }
        cerr<<"Too many steps in routine, no convergence\n";
        return ss;
    }

    //Utility routine used by vegas, to rebin a vector of densities xi into new bins defined by a
    //vector r.
    void rebin(Double_t rc, int nd, Double_t r[], Double_t xin[], Double_t xi[])
    {
        int i,k=0;
        Double_t dr=0.0,xn=0.0,xo=0.0;
        for (i=1;i<nd;i++) {
            while (rc > dr) dr += r[++k];
            if (k > 1) xo=xi[k-1];
            xn=xi[k];
            dr -= rc;
            xin[i]=xn-(xn-xo)*dr/r[k];
        }
        for (i=1;i<nd;i++) xi[i]=xin[i];
        xi[nd]=1.0;
    }


    //Performs Monte Carlo integration of a user-supplied ndim-dimensional function fxn over a
    //rectangular volume specified by regn[1..2*ndim], a vector consisting of ndim \lower left"
    //coordinates of the region followed by ndim \upper right" coordinates. The integration consists
    //of itmx iterations, each with approximately ncall calls to the function. After each iteration
    //the grid is refined; more than 5 or 10 iterations are rarely useful. The input flag init signals
    //whether this call is a new start, or a subsequent call for additional iterations (see comments
    //below). The input flag nprn (normally 0) controls the amount of diagnostic output. Returned
    //answers are tgral (the best estimate of the integral), sd (its standard deviation), and chi2a
    //(chisq2 per degree of freedom, an indicator of whether consistent results are being obtained). See
    //NR for further details.
    //init is whether initial or not
    //ncall is the number of function call
    //itmax is the max total number of interations done.
    //idum is used for random number initialization.
    //nprn is diagnostic flag turned off for the moment.
    //NOTE IS AS YET UNTESTED.
    void vegas(math_multidim_function *fxn, Double_t regn[], int ndim, int init,
    unsigned long ncall, int itmx, int nprn, long int idum, Double_t *tgral, Double_t *sd,
    Double_t *chi2a)
    {
    #define ALPH 1.5
    #define NDMX 50
    #define MXDIM 10
    #define TINY 1.0e-30
        //Best make everything static, allowing restarts.
        static int i,it,j,k,mds,nd,ndo,ng,npg,ia[MXDIM+1],kg[MXDIM+1];
        static Double_t calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo;
        static Double_t d[NDMX+1][MXDIM+1],di[NDMX+1][MXDIM+1],dt[MXDIM+1],
        dx[MXDIM+1], r[NDMX+1],x[MXDIM+1],xi[MXDIM+1][NDMX+1],xin[NDMX+1];
        static Double_t schi,si,swgt;
        //Normal entry. Enter here on a cold start.
        if (init <= 0) {
            mds=ndo=1; //Change to mds=0 to disable stratified sampling,
            for (j=1;j<=ndim;j++) xi[j][1]=1.0; //i.e., use importance sampling only.
        }
        //Enter here to inherit the grid from a previous call, but not its answers.
        if (init <= 1) si=swgt=schi=0.0;
        if (init <= 2) {// Enter here to inherit the previous grid and its answers
            nd=NDMX;
            ng=1;
            if (mds) { //Set up for stratification.
                ng=(int)pow(ncall/2.0+0.25,1.0/ndim);
                mds=1;
                if ((2*ng-NDMX) >= 0) {
                    mds = -1;
                    npg=ng/NDMX+1;
                    nd=ng/npg;
                    ng=npg*nd;
                }
            }
            for (k=1,i=1;i<=ndim;i++) k *= ng;
            npg=max(ncall/k,(long unsigned int)2);
            calls=(Double_t)npg * (Double_t)k;
            dxg=1.0/ng;
            for (dv2g=1,i=1;i<=ndim;i++) dv2g *= dxg;
            dv2g=pow(calls*dv2g,(Double_t)2.0)/npg/npg/(npg-1.0);
            xnd=nd;
            dxg *= xnd;
            xjac=1.0/calls;
            for (j=1;j<=ndim;j++) {
                dx[j]=regn[j+ndim]-regn[j];
                xjac *= dx[j];
            }
            if (nd != ndo) {//Do binning if necessary.
                for (i=1;i<=max(nd,ndo);i++) r[i]=1.0;
                for (j=1;j<=ndim;j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
                ndo=nd;
            }
			//comment this chunk
            /*if (nprn >= 0) {
                printf("%s: ndim= %3d ncall= %8.0f\n",
                " Input parameters for vegas",ndim,calls);
                printf("%28s it=%5d itmx=%5d\n"," ",it,itmx);
                printf("%28s nprn=%3d ALPH=%5.2f\n"," ",nprn,ALPH);
                printf("%28s mds=%3d nd=%4d\n"," ",mds,nd);
                for (j=1;j<=ndim;j++) {
                    printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
                    " ",j,regn[j],j,regn[j+ndim]);
                }
            }*/
        }
        for (it=1;it<=itmx;it++) {
            //Main iteration loop. Can enter here (init >= 3) to do an additional itmx iterations with
            //all other parameters unchanged.
            ti=tsi=0.0;
            for (j=1;j<=ndim;j++) {
                kg[j]=1;
                for (i=1;i<=nd;i++) d[i][j]=di[i][j]=0.0;
            }
            for (;;) {
                fb=f2b=0.0;
                for (k=1;k<=npg;k++) {
                    wgt=xjac;
                    for (j=1;j<=ndim;j++) {
                        xn=(kg[j]-ran2(&idum))*dxg+1.0;
#ifdef QUADPRECISION
                        ia[j]=max(min(to_int(xn),NDMX),1);
#elif QUADQUADPRECISION
                        ia[j]=max(min(to_int(xn),NDMX),1);
#else
                        ia[j]=max(min((int)(xn),NDMX),1);
#endif
                        if (ia[j] > 1) {
                            xo=xi[j][ia[j]]-xi[j][ia[j]-1];
                            rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
                        } else {
                            xo=xi[j][ia[j]];
                            rc=(xn-ia[j])*xo;
                        }
                        x[j]=regn[j]+rc*dx[j];
                        wgt *= xo*xnd;
                    }
                    //????? NR code as where wgt is weight
                    //f=wgt*(*fxn)(x,wgt);
                    //but! if f is an n dimensional function needs no weight
                    f=wgt*(fxn->function(x,ndim,fxn->params));
                    f2=f*f;
                    fb += f;
                    f2b += f2;
                    for (j=1;j<=ndim;j++) {
                        di[ia[j]][j] += f;
                        if (mds >= 0) d[ia[j]][j] += f2;
                    }
                }
                f2b=sqrt(f2b*npg);
                f2b=(f2b-fb)*(f2b+fb);
                if (f2b <= 0.0) f2b=TINY;
                ti += fb;
                tsi += f2b;
                if (mds < 0) {// Use stratified sampling.
                    for (j=1;j<=ndim;j++) d[ia[j]][j] += f2b;
                }
                for (k=ndim;k>=1;k--) {
                    kg[k] %= ng;
                    if (++kg[k] != 1) break;
                }
                if (k < 1) break;
            }
            tsi *= dv2g; //Compute final results for this iteration.
            wgt=1.0/tsi;
            si += wgt*ti;
            schi += wgt*ti*ti;
            swgt += wgt;
            *tgral=si/swgt;
            *chi2a=(schi-si*(*tgral))/(it-0.9999);
            if (*chi2a < 0.0) *chi2a = 0.0;
            *sd=sqrt(1.0/swgt);
            tsi=sqrt(tsi);
			//comment this chunk
            /*if (nprn >= 0) {
                printf("%s %3d : integral = %14.7g +/- %9.2g\n",
                " iteration no.",it,ti,tsi);
                printf("%s integral =%14.7g+/-%9.2g chi**2/IT n = %9.2g\n",
                " all iterations: ",*tgral,*sd,*chi2a);
                if (nprn) {
                    for (j=1;j<=ndim;j++) {
                        printf(" DATA FOR axis %2d\n",j);
                        printf("%6s%13s%11s%13s%11s%13s\n",
                        "X","delta i","X","delta i","X","delta i");
                        for (i=1+nprn/2;i<=nd;i += nprn+2) {
                            printf("%8.5f%12.4g%12.5f%12.4g%12.5f%12.4g\n",
                            xi[j][i],di[i][j],xi[j][i+1],
                            di[i+1][j],xi[j][i+2],di[i+2][j]);
                        }
                    }
                }
            }*/
            for (j=1;j<=ndim;j++) {
    // Refine the grid. Consult references to understand the subtlety of this procedure. 
    //The refinement is damped, to avoid rapid, destabilizing changes, and also compressed in range
    //by the exponent ALPH.
                xo=d[1][j];
                xn=d[2][j];
                d[1][j]=(xo+xn)/2.0;
                dt[j]=d[1][j];
                for (i=2;i<nd;i++) {
                    rc=xo+xn;
                    xo=xn;
                    xn=d[i+1][j];
                    d[i][j] = (rc+xn)/3.0;
                    dt[j] += d[i][j];
                }
                d[nd][j]=(xo+xn)/2.0;
                dt[j] += d[nd][j];
            }
            for (j=1;j<=ndim;j++) {
                rc=0.0;
                for (i=1;i<=nd;i++) {
                    if (d[i][j] < TINY) d[i][j]=TINY;
                    r[i]=pow((Double_t)((1.0-d[i][j]/dt[j])/(log(dt[j])-log(d[i][j]))),(Double_t)ALPH);
                    rc += r[i];
                }
                rebin(rc/xnd,nd,r,xin,xi[j]);
            }
        }
    }

    //integrate using vegas integrator.
    Double_t IntegrateVegasMonte(math_multidim_function * fm, Double_t *a,
                         Double_t *b, const int numintervals, Double_t *error, Double_t*chisq, int interations)
    {
        Double_t result,temperror,tempchisq;
        //construct regions array
        Double_t region[fm->ndim*2+1];
        for (int i=0;i<fm->ndim;i++){
            region[2*i+1]=a[i];
            region[2*i+1+1]=b[i];
        }
        vegas(fm, region, fm->ndim, 0, numintervals, interations, 0, lrand48(), &result, &temperror,&tempchisq);
        if (error)*error=temperror;
        if (chisq)*chisq=tempchisq;
        return result;
    }

    double IntegrateVegasMonte(gsl_monte_function * gslfm, double *a,
                         double *b, const int numintervals, double *error, double *chisq)
    {
        double  result,temperror;
        gsl_monte_vegas_state *ws = gsl_monte_vegas_alloc (gslfm->dim);
        const gsl_rng_type *T;
        gsl_rng *rnd;
        gsl_rng_env_setup ();
        T = gsl_rng_default;
        rnd = gsl_rng_alloc (T);
        gsl_monte_vegas_integrate (gslfm, a, b, gslfm->dim, numintervals, rnd, ws,&result, &temperror);
        if (error)*error=temperror;
        if (chisq)*chisq=ws->chisq;
        gsl_rng_free (rnd);
        gsl_monte_vegas_free (ws);
        return result;
    }

    //uses vegas integrator, sampling increased linearly until max (sampling*(n+1)) is reached.
    Double_t IntegrateRombergMonte(math_multidim_function * fm, Double_t *a,
                         Double_t *b, const Double_t epsrel, const int sampling, const int n, const int nextrap, int interations)
    {
        Double_t ss,dss;
        //stores trapezoidal approximations and stepsizes
        Double_t *s=new Double_t[n+1];
        Double_t *h=new Double_t[n+2];

        Double_t region[fm->ndim*2+1];
        for (int i=0;i<fm->ndim;i++){
            region[2*i+1]=a[i];
            region[2*i+1+1]=b[i];
        }

        //note that 
        h[0]=1.0;
        for (int j=0;j<n;j++)
        {
            int numintervals=sampling*(j+1);
            Double_t result, error, chisq;
            vegas(fm, region, fm->ndim, 0, numintervals, interations, 0, lrand48(), &result, &error,&chisq);
            /*gsl_monte_vegas_state *ws = gsl_monte_vegas_alloc (fm->dim);
            const gsl_rng_type *T;
            gsl_rng *rnd;
            gsl_rng_env_setup ();
            T = gsl_rng_default;
            rnd = gsl_rng_alloc (T);
            gsl_monte_vegas_integrate (fm, a, b, fm->dim, numintervals, rnd, ws,&result, &error);
            gsl_rng_free (rnd);
            gsl_monte_vegas_free (ws);*/
            s[j]=result;
            if (j>=nextrap-1){
                    PolyInt(&h[j-nextrap+1], &s[j-nextrap+1], nextrap, 0.0, ss, dss);
                    if (fabs(dss)<=epsrel*fabs(ss))return ss;
            }
            h[j+1]=0.25*h[j];
            //factor of 0.25 even though stepsize only decreaseed by 0.5 makes polynomial in extrapolate decrease in step size 
        }
        cerr<<"Too many steps in routine, no convergence\n";
        return ss;
    }

}
