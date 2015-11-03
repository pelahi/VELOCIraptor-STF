/*! \file Matrix2D.h
 *  \brief this file defines the \ref Math::Matrix2D class that is part of the Math namespace
 */

#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <iostream>
#include <complex>
#include <Precision.h>
#include <Coordinate2D.h>
#include <cmath>

using namespace std;
namespace Math
{
/*! \class Math::Matrix2D 
    \brief Basic 2x2 matrix.  
    Optimized implementation of 2x2 real matrix. Will eventually phase out and replace with more general matrix class \ref Math::GMatrix
*/
class Matrix2D
{
    private:
            
    Double_t matrix2d[4];
    
    public:
    
    /// \name Constructors
    //@{
    /// Constructor
    Matrix2D() { }
    
    /// Constructor: assign one number to all the elements
    Matrix2D(Double_t a)
    {
        for (int i = 0; i < 4; i++) matrix2d[i] = a;
    }
    
    /// Constructor: pass all 4 elements
    Matrix2D(Double_t a, Double_t b,
             Double_t c, Double_t d)
    {
        matrix2d[0] = a; matrix2d[1] = b;
        matrix2d[2] = c; matrix2d[3] = d;
    }
    
    /// Constructor: an array
    Matrix2D(Double_t m[4])
    {
        for (int i = 0; i < 4; i++) matrix2d[i] = m[i];
    }
    //@}

    /// \name Overloaded operators
    //@{
    /// How to access a matrix element: M(i, j)
    Double_t& operator () (int i, int j) { return matrix2d[j + 2 * i]; }
    const Double_t& operator() (int i, int j) const
    {
        return matrix2d[j + 2 * i];
    } 

    /// Matrix Addition
    Matrix2D operator + (const Matrix2D& m)
    {
        return Matrix2D(matrix2d[0] + m.matrix2d[0], 
                        matrix2d[1] + m.matrix2d[1],
                        matrix2d[2] + m.matrix2d[2], 
                        matrix2d[3] + m.matrix2d[3]);
    }
    
    /// Matrix Subtraction
    Matrix2D operator - (const Matrix2D& m)
    {
        return Matrix2D(matrix2d[0] - m.matrix2d[0], 
                        matrix2d[1] - m.matrix2d[1],
                        matrix2d[2] - m.matrix2d[2], 
                        matrix2d[3] - m.matrix2d[3]);
    }
    
    /// Multiply by a scalar: a * M   
    friend Matrix2D operator * (Double_t a, const Matrix2D& m);

    /// Multiply by a scalar: M * a
    Matrix2D operator * (Double_t a)
    {
        return Matrix2D(matrix2d[0] * a, matrix2d[1] * a,
                        matrix2d[2] * a, matrix2d[3] * a);
    }
    
    /// Multiply by a vector (2D coordinate): applicable
    /// under the new two-dimensional coordinate paradigm
    Coordinate2D operator * (const Coordinate2D& c)
    {
        return Coordinate2D(c.X() * matrix2d[0] + c.Y() * matrix2d[1],
                            c.X() * matrix2d[2] + c.Y() * matrix2d[3]);
    }                     
    
    /// Multiply by another matrix
    Matrix2D operator * (const Matrix2D& m)
    {
        return Matrix2D(matrix2d[0] * m.matrix2d[0] + 
                        matrix2d[1] * m.matrix2d[2],
                        matrix2d[0] * m.matrix2d[1] + 
                        matrix2d[1] * m.matrix2d[3],
                        matrix2d[2] * m.matrix2d[0] + 
                        matrix2d[3] * m.matrix2d[2],
                        matrix2d[2] * m.matrix2d[1] + 
                        matrix2d[3] * m.matrix2d[3]);
    }

    /// Print out the matrix
    friend std::ostream &operator << (std::ostream& stream,
                                      const Matrix2D& m);
    //@}
    
    /// \name Math operations
    //@{
    /// Transpose the matrix
    Matrix2D Transpose()
    {
        return Matrix2D(matrix2d[0], matrix2d[2],
                        matrix2d[1], matrix2d[3]);
    }
    
    /// Transpose matrix in place (changes the matrix permanently!)
    void TransposeInPlace()
    {
        Double_t temp;
        temp = matrix2d[1]; matrix2d[1] = matrix2d[2]; matrix2d[2] = temp;
    }
    
    /// Determinant of matrix
    Double_t Det()
    {
        return matrix2d[0] * matrix2d[3] - matrix2d[1] * matrix2d[2];
    }
    
    /// Adjugate (classical adjoint) of matrix
    Matrix2D Adjugate()
    {
        return Matrix2D( matrix2d[3], - matrix2d[1],
                        - matrix2d[2],  matrix2d[0]);
    }
    
    /// Inverse of matrix
    Matrix2D Inverse()
    {
        Double_t d = Det();
        if (d == 0.0)
            return Matrix2D(0, 0, 0, 0);
        else
            return Adjugate() * (1.0 / d);
    }
    
    /// Calculate the eigenvalues. Returns a two dimensional coordinate.
    Coordinate2D Eigenvalues() const;
    
    /// Compute the eigenvectors of the matrix if symmetric
    Matrix2D Eigenvectors(const Coordinate2D& e) const;
    //@}
    
};
}

#endif
