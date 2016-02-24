/*! \file Matrix.cxx
 *  \brief subroutines for \ref Math::Matrix class
 */

#include <Coordinate.h>
#include <Matrix.h>

using namespace std;
namespace Math
{
Matrix operator * (Double_t a, const Matrix& m)
{
    return Matrix(m.matrix[0] * a, m.matrix[1] * a, m.matrix[2] * a,
                  m.matrix[3] * a, m.matrix[4] * a, m.matrix[5] * a,
                  m.matrix[6] * a, m.matrix[7] * a, m.matrix[8] * a);
}


std::ostream &operator << (std::ostream& stream, const Matrix& m)
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            stream.width(6); stream.precision(3);
            stream << m.matrix[j + 3 * i] << " ";
        }
        stream << "\n";
    }
    return stream;
}
}
