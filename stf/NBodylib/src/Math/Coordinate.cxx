/*! \file Coordinate.cxx
 *  \brief subroutines of the \ref Math::Coordinate class that is part of the Math namespace
 */

#include <Coordinate.h>

using namespace std;
namespace Math
{
    Coordinate operator * (Double_t a, const Coordinate& c)
    {
        return Coordinate(a * c.coord[0], a * c.coord[1], a * c.coord[2]);
    }

    std::ostream& operator << (std::ostream &outs, Coordinate c)
    {
        outs.setf(std::ios::scientific);
        outs.width(16); outs.precision(7);
        outs << c.coord[0];
        outs.width(16); outs.precision(7);
        outs << c.coord[1];
        outs.width(16); outs.precision(7);
        outs << c.coord[2];
        return outs;
    }

    std::istream& operator >> (std::istream &ins, Coordinate& c)
    {   
        ins >> c.coord[0] >> c.coord[1] >> c.coord[2];
        return ins;
    }
}
