/*! \file SpecialFunctions.cxx
 *  \brief subroutines of special mathematical functions
 */

#include <SpecialFunctions.h>

using namespace std;
namespace Math
{
    int Factorial(int k)
    {
        int temp=1;
        for (int i=1;i<=k;i++) temp*=i;
        return temp;
    }

    //calculates spherical harmonics
/*    Double_t SphericalHarmonicR(Double_t theta, Double_t phi, int l, int m)
    {
        Double_t c=sqrt((Double_t)(2.0*l+1.0)/(Double_t)(4.0*M_PI)*(Double_t)Factorial(l-m)/(Double_t)Factorial(l+m));
        Double_t leg=AssociatedLegendrePolynomial();
    }
    Double_t SphericalHarmonicI(Double_t theta, Double_t phi, int l, int m)
    {
        Double_t c=sqrt((Double_t)(2.0*l+1.0)/(Double_t)(4.0*M_PI)*(Double_t)Factorial(l-m)/(Double_t)Factorial(l+m));
    }

    Double_t SphericalHarmonic(Double_t theta, Double_t phi, int l, int m)
    {
    }
*/
}
