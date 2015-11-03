/*! \file Coordinate2D.cxx
 *  \brief subroutines of the \ref Math::Coordinate2D class that is part of the Math namespace
 */

#include <Coordinate2D.h>

using namespace std;
namespace Math
{
    Coordinate2D operator * (Double_t a, const Coordinate2D& c)
    {
        return Coordinate2D(a * c.coord[0], a * c.coord[1]);
    }

    std::ostream& operator << (std::ostream &outs, Coordinate2D c)
    {
        outs.setf(std::ios::scientific);
        outs.width(16); outs.precision(7);
        outs << c.coord[0];
        outs.width(16); outs.precision(7);
        outs << c.coord[1];
        return outs;
    }

    std::istream& operator >> (std::istream &ins, Coordinate2D& c)
    {
        ins >> c.coord[0] >> c.coord[1];
        return ins;
    }
}
