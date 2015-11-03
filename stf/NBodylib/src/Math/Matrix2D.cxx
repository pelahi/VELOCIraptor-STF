/*! \file Matrix2D.cxx
 *  \brief subroutines for \ref Math::Matrix2D class
 */

#include <Coordinate2D.h>
#include <Matrix2D.h>

using namespace std;
namespace Math
{
Matrix2D operator * (Double_t a, const Matrix2D& m)
{
    return Matrix2D(m.matrix2d[0] * a, m.matrix2d[1] * a, 
                    m.matrix2d[2] * a, m.matrix2d[3] * a);
}


std::ostream &operator << (std::ostream& stream, const Matrix2D& m)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            stream.width(6); stream.precision(3);
            stream << m.matrix2d[j + 2 * i] << " ";
        }
        if (i != 2) stream << "\n";
    }
    return stream;
}

// Calculate the eigenvalues of a two-dimensional matrix. Stored
// as a two dimensional vector e1, e2. The eigenvalues may be real
// or complex; in the complex case the eigenvalues are conjugate
// pairs, in which case the returned coordinate will consist of
// the real part for e1 and the imaginary part for e2.
Coordinate2D Matrix2D::Eigenvalues() const
{
    // First of all, if the matrix is symmetric, the 
    // eigenvalues are real, as they should be. If the matrix is
    // not symmetric, things are slightly more complicated.
    // The discriminant of the characteristic polynomial
    Double_t D = (matrix2d[0] + matrix2d[3]) * (matrix2d[0] + matrix2d[3]) -
               4.0 * matrix2d[0] * matrix2d[3] +
               4.0 * matrix2d[1] * matrix2d[2];
    if (D > 0.0)
    {
        Double_t e1 = 0.5 * (matrix2d[0] + matrix2d[3] + sqrt(D));
        Double_t e2 = 0.5 * (matrix2d[0] + matrix2d[3] - sqrt(D));
        
        return Coordinate2D(e1, e2);
    }
    // Else this function returns the real and imaginary parts separately
    else // (D < 0.0)
    {
        Double_t e1 = 0.5 * (matrix2d[0] + matrix2d[3]);
        Double_t e2 = 0.5 * (sqrt(-D));
        
        return Coordinate2D(e1, e2);
    }
}

// Compute the eigenvectors of the matrix. The following code 
// works for symmetric 2 x 2 matrices provided that neither off-diagonal
// element is zero and assumes a non-singular matrix
Matrix2D Matrix2D::Eigenvectors(const Coordinate2D& e) const
{
    Double_t g[2];
    Matrix2D m;
    
    if (matrix2d[1] == 0 && matrix2d[2] == 0)
    {
        if (matrix2d[0] > matrix2d[3])
        {
            m(0, 0) = m(1, 1) = 1.0;
            m(1, 0) = m(0, 1) = 0.0;
        }
        else
        {
            m(0, 0) = m(1, 1) = 0.0;
            m(1, 0) = m(0, 1) = 1.0;
        }
    }
    else if (matrix2d[1] == 0 && matrix2d[2] != 0)
    {
        if (matrix2d[0] > matrix2d[3])
        {
            Double_t h = (matrix2d[3] - matrix2d[0]) / matrix2d[2];
            m(0, 0) = - h / sqrt(h * h + 1.0);
            m(0, 1) = 1.0 / sqrt(h * h + 1.0);
            m(1, 0) = 0.0;
            m(1, 1) - 1.0;
        }
        else
        {
            Double_t h = (matrix2d[3] - matrix2d[0]) / matrix2d[2];
            m(1, 0) = - h / sqrt(h * h + 1.0);
            m(1, 1) = 1.0 / sqrt(h * h + 1.0);
            m(0, 0) = 0.0;
            m(0, 1) = 1.0;               
        }
    }
    else if (matrix2d[1] != 0 && matrix2d[2] == 0)
    {
        if (matrix2d[0] > matrix2d[3])
        {
            Double_t h =  matrix2d[1] / (matrix2d[0] - matrix2d[3]);
            m(1, 0) = - h / sqrt(h * h + 1.0);
            m(1, 1) = 1.0 / sqrt(h * h + 1.0);
            m(0, 1) = 0.0;
            m(0, 0) = 1.0;
        }
        else
        {
            Double_t h = matrix2d[1] / (matrix2d[0] - matrix2d[3]);
            m(0, 0) = - h / sqrt(h * h + 1.0);
            m(0, 1) = 1.0 / sqrt(h * h + 1.0);
            m(1, 1) = 0.0;
            m(1, 0) = 1.0;               
        }
    }
    else
    {
        for (int i = 0; i < 2; i++)
        {
            g[i] = - matrix2d[1] / (matrix2d[0] - e[i]);
            m(i, 0) = g[i] / sqrt(g[i] * g[i] + 1.0);
            m(i, 1) = 1.0 / sqrt(g[i] * g[i] + 1.0);
        }
    }
    
    return m;
}
}
