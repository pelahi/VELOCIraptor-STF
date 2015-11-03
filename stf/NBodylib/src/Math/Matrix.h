/*! \file Matrix.h
 *  \brief this file defines the \ref Math::Matrix class that is part of the Math namespace
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <complex>
#include <Precision.h>
#include <Coordinate.h>
#include <cmath>
#include <cstdio>

using namespace std;
namespace Math
{

/*! \class Math::Matrix 
    \brief Basic 3x3 matrix.  
    Optimized implementation of 3x3 real matrix. Will eventually phase out and replace with more general matrix class \ref Math::GMatrix
*/
class Matrix
{
    private:

    Double_t matrix[9];

    public:

    ///\name Constructors
    //@{
    ///Constructor (empty).
    Matrix() { }
    /// Constructor: assign one number to all elements.
    Matrix(Double_t a)
    {
        for (int i = 0; i < 9; i++) matrix[i] = a;
    }
    /// Constructor: pass all 9 elements.
    Matrix(Double_t a, Double_t b, Double_t c,
           Double_t d, Double_t e, Double_t f,
           Double_t g, Double_t h, Double_t i)
    {
        matrix[0] = a; matrix[1] = b; matrix[2] = c;
        matrix[3] = d; matrix[4] = e; matrix[5] = f;
        matrix[6] = g; matrix[7] = h; matrix[8] = i;
    }
    /// Constructor: an array (in row-major format).
    Matrix(const Double_t m[9])
    {
        for (int i = 0; i < 9; i++) matrix[i] = m[i];
    }
    /// Copy Constructor
    Matrix(const Matrix &m)
    {
        for (int i = 0; i < 9; i++) matrix[i] = m.matrix[i];
    }
    //@}
    
    /// \name Overloaded operators
    //@{
    /// How to access a matrix element: M(i, j)
    Double_t& operator () (int i, int j) { return matrix[j + 3 * i]; }
    const Double_t& operator() (int i, int j) const
    {
        return matrix[j + 3 * i];
    }

    /// Matrix addition
    inline Matrix operator + (const Matrix& m)
    {
        return Matrix(matrix[0] + m.matrix[0], matrix[1] + m.matrix[1],
                      matrix[2] + m.matrix[2], matrix[3] + m.matrix[3],
                      matrix[4] + m.matrix[4], matrix[5] + m.matrix[5],
                      matrix[6] + m.matrix[6], matrix[7] + m.matrix[7],
                      matrix[8] + m.matrix[8]);
    }

    /// Matrix subtraction
    inline Matrix operator - (const Matrix& m)
    {
        return Matrix(matrix[0] - m.matrix[0], matrix[1] - m.matrix[1],
                      matrix[2] - m.matrix[2], matrix[3] - m.matrix[3],
                      matrix[4] - m.matrix[4], matrix[5] - m.matrix[5],
                      matrix[6] - m.matrix[6], matrix[7] - m.matrix[7],
                      matrix[8] - m.matrix[8]);
    }

    /// Multiply by a scalar: a * M
    friend Matrix operator * (Double_t a, const Matrix& m);

    /// Multiply by a scalar: M * a
    inline Matrix operator * (Double_t a)
    {
        return Matrix(matrix[0] * a, matrix[1] * a, matrix[2] * a,
                      matrix[3] * a, matrix[4] * a, matrix[5] * a,
                      matrix[6] * a, matrix[7] * a, matrix[8] * a);
    }

    /// Multiply by a vector (coordinate)
    inline Coordinate operator * (const Coordinate& c)
    {
        return Coordinate(c.X() * matrix[0] + c.Y() * matrix[1] + 
                          c.Z() * matrix[2],  c.X() * matrix[3] + 
                          c.Y() * matrix[4] + c.Z() * matrix[5],
                          c.X() * matrix[6] + c.Y() * matrix[7] + 
                          c.Z() * matrix[8]);
    }

    /// Multiply by another matrix
    inline Matrix operator * (const Matrix& m)
    {
        return Matrix(matrix[0] * m.matrix[0] + matrix[1] * m.matrix[3] +
                      matrix[2] * m.matrix[6],
                      matrix[0] * m.matrix[1] + matrix[1] * m.matrix[4] +
                      matrix[2] * m.matrix[7],
                      matrix[0] * m.matrix[2] + matrix[1] * m.matrix[5] +
                      matrix[2] * m.matrix[8],
                      matrix[3] * m.matrix[0] + matrix[4] * m.matrix[3] +
                      matrix[5] * m.matrix[6],
                      matrix[3] * m.matrix[1] + matrix[4] * m.matrix[4] +
                      matrix[5] * m.matrix[7],
                      matrix[3] * m.matrix[2] + matrix[4] * m.matrix[5] +
                      matrix[5] * m.matrix[8],
                      matrix[6] * m.matrix[0] + matrix[7] * m.matrix[3] +
                      matrix[8] * m.matrix[6],
                      matrix[6] * m.matrix[1] + matrix[7] * m.matrix[4] +
                      matrix[8] * m.matrix[7],
                      matrix[6] * m.matrix[2] + matrix[7] * m.matrix[5] +
                      matrix[8] * m.matrix[8]);
    }

    /// Asignment constructor
    Matrix& operator = (const Matrix m)
    {
        for (int i = 0; i < 9; i++) matrix[i] = m.matrix[i];
        return *this;
    }
    //@}
    
    /// \name Math operations 
    //@{
    /// Trace the matrix.
    Double_t Trace() const
    {
        return matrix[0]*matrix[4]*matrix[8];
    }

    /// Transpose the matrix.
    Matrix Transpose() const
    {
        return Matrix(matrix[0], matrix[3], matrix[6],
                      matrix[1], matrix[4], matrix[7],
                      matrix[2], matrix[5], matrix[8]);
    }

    /// Transope matrix in place (changes the matrix permanently!)
    void TransposeInPlace()
    {
        Double_t temp;
        temp = matrix[1];  matrix[1] = matrix[3]; matrix[3] = temp;
        temp = matrix[2];  matrix[2] = matrix[6]; matrix[6] = temp;
        temp = matrix[5];  matrix[5] = matrix[7]; matrix[7] = temp;
    }

    /// Determinant of matrix
    inline Double_t Det() const
    {
        return matrix[0] * (matrix[4] * matrix[8] - matrix[5] * matrix[7])
             - matrix[1] * (matrix[3] * matrix[8] - matrix[5] * matrix[6])
             + matrix[2] * (matrix[3] * matrix[7] - matrix[4] * matrix[6]);
    }

    /// Adjugate (classical adjoint) of matrix
    inline Matrix Adjugate() const 
    {
        return Matrix( (matrix[4] * matrix[8] - matrix[5] * matrix[7]),
                       (matrix[2] * matrix[7] - matrix[1] * matrix[8]),
                       (matrix[1] * matrix[5] - matrix[2] * matrix[4]),
                       (matrix[5] * matrix[6] - matrix[3] * matrix[8]),
                       (matrix[0] * matrix[8] - matrix[2] * matrix[6]),
                       (matrix[2] * matrix[3] - matrix[0] * matrix[5]),
                       (matrix[3] * matrix[7] - matrix[4] * matrix[6]),
                       (matrix[1] * matrix[6] - matrix[0] * matrix[7]),
                       (matrix[0] * matrix[4] - matrix[1] * matrix[3]));
    }

    /*! \brief Lower-Upper decompistion
    
        LU decomposition of matrix stored in such a fashion that 
        matrix return is (for 3x3): 
        \f[
        \left(\begin{array}{ccc}
            L11 & U12 &U13\\
            L21 & L22 & U23\\
            L31 & L32 & L33
        \end{array}\right)
        \f]
        L is lower triangular matrix, U is upper unit triangular matrix. \n
        Crout's algorithm is used to determine the L and U matrices. System of equations to solve are : \n 
        \f[ 
        \begin{array}{lcl}
            i<j : & l_{i1}*u_{1j}+l_{i2}*u_{2j}+...+l_{ii}*u_{ij}&=a_{ij} \\
            i=j : & l_{i1}*u_{1j}+l_{i2}*u_{2j}+...+l_{ii}*u_{jj}&=a_{ij} \\
            i>j : & l_{i1}*u_{1j}+l_{i2}*u_{2j}+...+l_{ij}*u_{ij}&=a_{ij} 
            \end{array}
        \f]
        with condition that \f$ l_{ii}=1 \f$. \n 
        If matrix singular and no decomposition, print error and return zero matrix. 
        Here matrix is singular if diagonal element is zero. \n  
        Must alter this to pivot method so that subroutine reorder matrix if necessary and diagonals
        can be zero. 
    */
    inline Matrix LUDecomp() const
    {
        Double_t A[9];
        int i, j, k, p, n=3;
        Double_t *p_k, *p_row, *p_col;
        for (j=0;j<9;j++)A[j]=matrix[j];
        for (k = 0, p_k = A; k < n; p_k += n, k++) {
            for (i = k, p_row = p_k; i < n; p_row += n, i++) {
                for (p = 0, p_col = A; p < k; p_col += n, p++)
                    *(p_row + k) -= *(p_row + p) * *(p_col + k);
            }
            if ( *(p_k + k) == 0.0 ) {
                printf("Singular matrix, returning zero matrix\n");
                return Matrix(0.);
            }
            for (j = k+1; j < n; j++) {
                for (p = 0, p_col = A; p < k; p_col += n,  p++)
                    *(p_k + j) -= *(p_k + p) * *(p_col + j);
                *(p_k + j) /= *(p_k + k);
            }
        }
        Matrix LU(A);
        return LU;
    }


    /// Inverse of matrix
    Matrix Inverse() const
    {
        Double_t d = Det();
        if (d==0.0) return Matrix(0, 0, 0, 0, 0, 0, 0, 0, 0);
        else return Adjugate() * (1.0 / d);
    }
    /// Print out the matrix
    friend std::ostream &operator << (std::ostream& stream, const Matrix& m);

    /*! \brief Calculate eigenvalues of matrix.  
    
        Eigenvalues stored as a Coordinate e1, e2, e3 ordered such that e1>e2>e3. 
        The user must know what to expect from this function, whether three real eigenvalues
        or one real and two complex (which are conjugate pairs).
        Only three numbers can be returned, so it's either all three
        real values, or one real value and two values that are the
        real and imaginary parts of the two complex values.
    */
    Coordinate Eigenvalues() const
    {
    // The eigenvalues y are the solutions to the characteristic 
    // equation det(M - yI) = 0.  det(M - yI) is a cubic equation.
    // First compute co-efficients of characteristic equation
    // y^3 + a_2 * y^2 + a_1 * y + a_0
    Double_t a2,a1,a0;
    //double a2,a1,a0;
    a2 = -matrix[0] - matrix[4] - matrix[8];
    a1 =  matrix[4] * matrix[8] - matrix[5] * matrix[7] +
                 matrix[0] * matrix[8] + matrix[0] * matrix[4] -
                 matrix[1] * matrix[3] - matrix[2] * matrix[6];
    a0 =  matrix[0] * matrix[5] * matrix[7] - 
                 matrix[0] * matrix[4] * matrix[8] +
                 matrix[1] * matrix[3] * matrix[8] -
                 matrix[1] * matrix[5] * matrix[6] -
                 matrix[2] * matrix[3] * matrix[7] +
                 matrix[2] * matrix[4] * matrix[6];

    // Following the solution given at
    // http://mathworld.wolfram.com/CubicEquation.html
    Double_t Q,R,D;
    //double Q,R,D;
    Q = (3.0 * a1 - a2 * a2) / 9.0;
    R = (9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2) / 54.0;

    // Some more convienience constants.
    D = Q * Q * Q + R * R;
    // Now, if D < 0, the eigenvalues will be all real and
    // different.  If D == 0, two eigenvalues will be the same.
    // Define unit imaginary number.

    std::complex<Double_t> I(0, 1);
    //std::complex<double> I(0, 1);
    if (D < 0.0)
    {
        std::complex<Double_t> S,T;
        //std::complex<double> S,T;

        S = pow(std::complex<Double_t>(R, sqrt(-D)),(Double_t)(1.0 / 3.0));
        T = pow(std::complex<Double_t>(R,-sqrt(-D)),(Double_t)(1.0 / 3.0));
        //S = pow(std::complex<double>(R, sqrt(-D)),(double)(1.0 / 3.0));
        //T = pow(std::complex<double>(R,-sqrt(-D)),(double)(1.0 / 3.0));

        // Here's the eigenvalues:
        Double_t e1,e2,e3;
        e1 = real(-(Double_t)(1.0/3.0)*a2 + S + T);
        e2 = real(-(Double_t)(1.0/3.0)*a2-(Double_t)0.5*(S + T)-(Double_t)0.5*I*(Double_t)sqrt(3.0)*(S - T));
        e3 = real(-(Double_t)(1.0/3.0)*a2-(Double_t)0.5*(S + T)+(Double_t)0.5*I*(Double_t)sqrt(3.0)*(S - T));
        //e2 = real(-1.0/3.0*a2-(double)0.5*(S + T)-(double)0.5*I*(double)sqrt(3.0)*(S - T));
        //e3 = real(-1.0/3.0*a2-(double)0.5*(S + T)+(double)0.5*I*(double)sqrt(3.0)*(S - T));

        return Coordinate(e1, e2, e3);
    }
    else if (D == 0.0)
    {
        Double_t S;

        // Handle cube root
        if (R < 0.0) S = -pow(-R, (Double_t)(1.0 / 3.0));
        else S = pow(R, (Double_t)(1.0 / 3.0));

        // eigenvalues, e2 and e3 the same
        Double_t e1,e2;
        e1 = -(Double_t)(1.0/3.0) * a2 + 2.0 * S;
        e2 = -(Double_t)(1.0/3.0) * a2 - S;

        return Coordinate(e1, e2, e2);
    }
    // If D > 0, one eigenvalue will be real, the other two will
    // be complex conjugate pairs.  That is, we'll have:
    // e1 = a,   e2 = b + ic,   e3 = b - ic.  Returned
    // is a, b, c.
    else // (D > 0.0)
    {
        // Have to deal with cube root properly.
        Double_t temp = R + sqrt(D);
        std::complex<Double_t> S;
        std::complex<Double_t> T;
        if (temp >= 0.0)
            S = pow(temp, (Double_t)(1.0 / 3.0));
        else
            S = -pow(-temp, (Double_t)(1.0 / 3.0));
        temp = R - sqrt(D);
        if (temp >= 0.0)
            T = pow(temp, (Double_t)(1.0 / 3.0));
        else
            T = -pow(-temp, (Double_t)(1.0 / 3.0));

        std::complex<Double_t> e1 = -(Double_t)(1.0/3.0) * a2 + S + T;
        std::complex<Double_t> e2 = -(Double_t)(1.0/3.0) * a2 - (Double_t)0.5 * (S + T) +
                                   (Double_t)0.5 * I * (Double_t)sqrt(3.0) * (S - T);

        return Coordinate(e1.real(), e2.real(), e2.imag());
    }
    }

    /*! \brief Computes the eigenvectors of the matrix, given the three eigenvalues.  
        
        Sadly, this function will only work properly on symmetric matrices -- i.e., all real eigenvalues.
        The matrix returned is 3 row vectors.
    */
    Matrix Eigenvectors(const Coordinate& e) const
    {
    Double_t g[3];
    Double_t h[3];
    Matrix m;
    for (int i = 0; i < 3; i++)
    {
        // If two eigenvalues are the same, deal with
        // differently (courtesy David Puglielli)
        if((e[i % 3] == e[(i + 1) % 3]) && (e[i % 3] != e[(i + 2) % 3]))
        {
            g[i] = matrix[1] / (matrix[0] - e[i]);
            h[i] = matrix[2] / (matrix[0] - e[i]);
            m(i % 3, 0) = - g[i] / sqrt(g[i] * g[i] + 1.0);
            m(i % 3, 1) = 1.0 / sqrt(g[i] * g[i] + 1.0);
            m(i % 3, 2) = 0.0;
            m((i + 1) % 3, 0) = - h[i] / sqrt(h[i] * h[i] + 1.0);
            m((i + 1) % 3, 1) = 0.0;
            m((i + 1) % 3, 2) = 1.0 / sqrt(h[i] * h[i] + 1.0);
        }
        else if (matrix[1] * matrix[1] - 
                (matrix[0] - e[i]) * (matrix[4] - e[i]) == 0.0)
        {
            // assumes that (matrix[5]*matrix[1] - matrix[2]*(matrix[4] - e[i]) == 0.0)
            g[i] = matrix[5] / matrix[2];
            m(i, 2) = 0.0;
            // Warning: assumes matrix[8]*matrix[3] - matrix[5]*matrix[6] != 0. I think this is
            // reasonable since then we would get two eigenvectors, but we already have code to 
            // account for that (although even that code is not completely general).
            m(i, 1) = 1.0 / sqrt(g[i] * g[i] + 1.0);
            m(i, 0) = - g[i] / sqrt(g[i] * g[i] + 1.0);
        }

        // Now deal with cases where some elements are zero
        else if(matrix[1] == 0.0 && matrix[2] == 0.0 && matrix[5] != 0.0 &&
                matrix[0] != 0.0 && matrix[4] != 0.0 && matrix[8] != 0.0)
        {
            g[i] = - matrix[5] / matrix[4];
            m(i, 1) = g[i] / sqrt(g[i] * g[i] + 1.0);
            m(i, 2) = 1.0 / sqrt(g[i] * g[i] + 1.0);
            m(i, 0) = 0.0;
        }
        else if(matrix[1] != 0.0 && matrix[2] == 0.0 && matrix[5] == 0.0 &&
                matrix[0] != 0.0 && matrix[4] != 0.0 && matrix[8] != 0.0)
        {
            g[i] = - matrix[1] / matrix[0];
            m(i, 0) = g[i] / sqrt(g[i] * g[i] + 1.0);
            m(i, 1) = 1.0 / sqrt(g[i] * g[i] + 1.0);
            m(i, 2) = 0.0;
        }
        else if(matrix[1] == 0.0 && matrix[2] != 0.0 && matrix[5] == 0.0 &&
                matrix[0] != 0.0 && matrix[4] != 0.0 && matrix[8] != 0.0)
        {
            g[i] = - matrix[2] / matrix[0];
            m(i, 0) = g[i] / sqrt(g[i] * g[i] + 1.0);
            m(i, 2) = 1.0 / sqrt(g[i] * g[i] + 1.0);
            m(i, 1) = 0.0;
        }
        else if(matrix[1] != 0.0 && matrix[2] != 0.0 && matrix[5] != 0.0 &&
                matrix[0] == 0.0 && matrix[4] == 0.0 && matrix[8] != 0.0)
        {
            g[i] = - matrix[5] / matrix[1];
            h[i] = - matrix[2] / matrix[1];
            m(i, 0) = g[i] / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
            m(i, 1) = h[i] / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
            m(i, 2) = 1.0 / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
        }
        else if(matrix[1] != 0.0 && matrix[2] != 0.0 && matrix[5] != 0.0 &&
                matrix[0] == 0.0 && matrix[4] != 0.0 && matrix[8] == 0.0)
        {
            g[i] = - matrix[5] / matrix[1];
            h[i] = - matrix[2] / matrix[1];
            m(i, 0) = - g[i] / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
            m(i, 1) = h[i] / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
            m(i, 2) = 1.0 / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
        }
        else if(matrix[1] != 0.0 && matrix[2] != 0.0 && matrix[5] != 0.0 &&
                matrix[0] != 0.0 && matrix[4] == 0.0 && matrix[8] == 0.0)
        {
            g[i] = - matrix[5] / matrix[2];
            h[i] = - matrix[5] / matrix[1];
            m(i, 1) = g[i] / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
            m(i, 2) = h[i] / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
            m(i, 0) = 1.0 / sqrt(g[i] * g[i] + h[i] * h[i] + 1.0);
        }
        else
        {
            g[i] = -(matrix[2] * (matrix[4] - e[i]) - matrix[1] * 
                     matrix[5]) / (matrix[1] * matrix[1] -
                    (matrix[0] - e[i]) * (matrix[4] - e[i]));
            h[i] = (matrix[1] * matrix[2] - (matrix[0] - e[i]) * 
                    matrix[5]) / (matrix[1] * matrix[1] - 
                   (matrix[0] - e[i]) * (matrix[4] - e[i]));
            m(i, 0) = - g[i] / sqrt (g[i] * g[i] + 
                                     h[i] * h[i] + 1.0);
            m(i, 1) = - h[i] / sqrt (g[i] * g[i] + 
                                     h[i] * h[i] + 1.0);
            m(i, 2) =    1.0 / sqrt (g[i] * g[i] + 
                                     h[i] * h[i] + 1.0);
        }
    }

    return m;
    }
    //@}
};
}
#endif
