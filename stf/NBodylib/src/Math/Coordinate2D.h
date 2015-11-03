/*! \file Coordinate2D.h
 *  \brief this file defines the \ref Math::Coordinate2D class that is part of the Math namespace
 */

#ifndef COORDINATE2D_H
#define COORDINATE2D_H

#include <iostream>
#include <cmath>
#include <Precision.h>

using namespace std;
namespace Math
{
/*!
    \class Math::Coordinate2D 
    \brief a simple 2-dimensional coordinate in real space.
*/
    class Coordinate2D
    {
        private:
    
        // the data: an array       
        Double_t coord[2];
    
        public:
    
        /// \name Constructors
        //@{
        ///from 2 numbers ...
        Coordinate2D(const Double_t x = 0, const Double_t y = 0)
        {
            coord[0] = x;
            coord[1] = y;
        }
        /// or from an array.  There is no bounds check here, so pass an
        /// array with at least 2 values.
        Coordinate2D(const Double_t *array)
        {
            coord[0] = array[0];
            coord[1] = array[1];
        }
        Coordinate2D(const Real_t *array)
        {
            coord[0] = array[0];
            coord[1] = array[1];
        }
        //@}

        /// \name Get routines
        //@{
        /// Look at x coordinate
        Double_t X() const { return coord[0]; }
        /// Look at y coordinate
        Double_t Y() const { return coord[1]; }
        /// Get at the array itself
        Double_t* Coord() { return coord; }
        //@}

        /// \name Overloaded operators
        /// \todo still to implement: +=, -=, *=, /=, ... ?
        //@{
        /// Assignment for an array.
        Coordinate2D& operator = (const Double_t *array)
        {
            coord[0] = array[0];
            coord[1] = array[1];
            return *this;
        }
        Coordinate2D& operator = (const Real_t *array)
        {
            coord[0] = array[0];
            coord[1] = array[1];
            return *this;
        }
        /// Interface like an array.  Can set values with this, too.
        Double_t& operator [] (int i) { return coord[i]; }
        /// Const version.
        const Double_t& operator [] (int i) const { return coord[i]; }
        /// Multiply by a scalar
        Coordinate2D operator * (Double_t a)
        {
            return Coordinate2D(a * coord[0], a * coord[1]);
        }
    
        /// Divide by a scalar
        Coordinate2D operator / (Double_t a)
        {
            return Coordinate2D(coord[0] / a, coord[1] / a);
        }
    
        /// Dot product -- coordinate times coordinate.
        Double_t operator * (const Coordinate2D& c)
        {
            return coord[0] * c.coord[0] + coord[1] * c.coord[1];
        }
    
        /// Adds two coordinates
        Coordinate2D operator + (const Coordinate2D& c)
        {
            return Coordinate2D(coord[0] + c.coord[0], coord[1] + c.coord[1]);
        }
    
        /// Subtracts two coordinates.
        Coordinate2D operator - (const Coordinate2D& c)
        {
            return Coordinate2D(coord[0] - c.coord[0], coord[1] - c.coord[1]);
        }
    
        /// Multiply by a scalar.
        friend Coordinate2D operator * (Double_t a, const Coordinate2D& c);

        /// Output to stream.
        friend std::ostream &operator << (std::ostream &outs, Coordinate2D c);
        /// Input from stream.
        friend std::istream &operator >> (std::istream &ins, Coordinate2D& c);

        //@}
        
        /// \name Math operations
        //@{
        /// Length of the vector (distance from origin to point).
        Double_t Length()
        {
            return sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
        }
    
        /// Normalized vector (unit vector along direction of vector).
        Coordinate2D Normal()
        {
            Double_t L = Length();
            return Coordinate2D(coord[0] / L, coord[1] / L);
        }
        //@}
    
    };
}

#endif    
