/*! \file Interpolate.cxx
 *  \brief interpolation subroutines 
 */

#include <Interpolate.h>

using namespace std;

namespace Math
{
    //polynomial interpolation
    /*void PolyInt(Double_t *xa, Double_t *ya, int n, Double_t x, Double_t &y, Double_t &dy)
    {
        const int nmax=10;
        int i,m,ns;
        Double_t den,dif,dift,ho,hp,w;
        Double_t *c=new Double_t[nmax],*d=new Double_t[nmax];

        ns=0;
        dif=fabs(x-xa[0]);
        for(int i=1;i<n;i++)
        {
            dift=fabs(x-xa[i]);
            if (dift<dif)
            {
                ns=i;
                dif=dift;
            }
            c[i]=ya[i];
            d[i]=ya[i];
        }
		y=ya[ns];

		cout<<"first estimate"<<y<<" "<<dy<<endl;
        for (int m=0;m<n-1;m++)
        {
            for (int i=1;i<=n-m;i++)
            for (int i=0;i<n-m;i++)
            {
                ho=xa[i]-x;
                hp=xa[i+m]-x;
                w=c[i+1]-d[i];
                den=ho-hp;
                if(den==0.){cerr<<"failure in polint\n";exit(8);}
                den=w/den;
                //d[i-1]=hp*den;
                //c[i-1]=ho*den;
                d[i]=hp*den;
                c[i]=ho*den;
            }
            if (2*ns<n-m) dy=c[ns+1];
            else
            {
                dy=d[ns--];
            }
            y=y+dy;
			//y+=(dy=(2*ns<(n-m)?c[ns+1]:d[ns--]));
        }
		cout<<"latestest estimate"<<y<<" "<<dy<<endl;
        delete []c;
        delete []d;
        return;
    }*/

    void PolyInt(Double_t *xa, Double_t *ya, int n, Double_t x, Double_t &y, Double_t &dy)
    {
        int ns;
        Double_t den,dif,dift,ho,hp,w;
        Double_t *c=new Double_t[n],*d=new Double_t[n];
        //interpolate uses index one arrays
        xa--;ya--;c--;d--;
        ns=1;
        dif=fabs(x-xa[1]);
        for(int i=1;i<=n;i++)
        {
            dift=fabs(x-xa[i]);
            if (dift<dif)
            {
                ns=i;
                dif=dift;
            }
            c[i]=ya[i];
            d[i]=ya[i];
        }
        y=ya[ns--];
        for (int m=1;m<n;m++)
        {
            for (int i=1;i<=n-m;i++)
            {
                ho=xa[i]-x;
                hp=xa[i+m]-x;
                w=c[i+1]-d[i];
                den=ho-hp;
                if(den==0.){cerr<<"failure in polint\n";exit(8);}
                den=w/den;
                d[i]=hp*den;
                c[i]=ho*den;
            }
            if (2*ns<n-m) dy=c[ns+1];
            else
            {
                dy=d[ns--];
            }
            y=y+dy;
        }
        //reset pointers
		c++;d++;
        delete []c;
        delete []d;
        xa++;ya++;
        return;
    }

}
