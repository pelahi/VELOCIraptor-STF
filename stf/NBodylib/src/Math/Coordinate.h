/*! \file Coordinate.h
 *  \brief this file defines the \ref Math::Coordinate class that is part of the Math namespace
 */


#ifndef COORDINATE_H
#define COORDINATE_H

#include <iostream>
#include <cmath>
#include <Precision.h>

using namespace std;
namespace Math
{

/// \struct Coord
/// \brief this structure very simple structure, 1 dim array with 3 cells. storage is small as essentially no functionality
struct Coord
{
    Double_t pos[3];

    //    OPERATORS
    Coord &operator=(const Coord &q)
    {
        pos[0]=q.pos[0];
        pos[1]=q.pos[1];
        pos[2]=q.pos[2];
        return *this;
    }

    Coord &operator=(const Double_t q[3])
    {
        pos[0]=q[0];
        pos[1]=q[1];
        pos[2]=q[2];
        return *this;
    }

};

/*!
    \class Math::Coordinate 
    \brief a simple 3-dimensional coordinate in real space.
*/
class Coordinate
{
    private:
    // the data: an array
    Double_t coord[3];

    public:

    /// \name Constructors
    //@{
    
    ///from 3 numbers ...
    Coordinate(const Double_t x = 0, const Double_t y = 0, const Double_t z = 0)
    {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    /// or from an array.  There is no bounds check here, so pass an array with at least 3 values.
    Coordinate(const Double_t *array)
    {
        coord[0] = array[0];
        coord[1] = array[1];
        coord[2] = array[2];
    }
    Coordinate(const Real_t *array)
    {
        coord[0] = array[0];
        coord[1] = array[1];
        coord[2] = array[2];
    }
    /// from \ref Math::Coord structure
    Coordinate(const Coord q)
    {
        coord[0] = q.pos[0];
        coord[1] = q.pos[1];
        coord[2] = q.pos[2];
    }
    /// copy constructor
    Coordinate(const Coordinate &q)
    {
        coord[0] = q[0];
        coord[1] = q[1];
        coord[2] = q[2];
    }

    //@}

    /// \name Get routines
    //@{
    /// Look at x coordinate
    Double_t X() const { return coord[0]; }
    /// Look at y coordinate
    Double_t Y() const { return coord[1]; }
    /// Look at z coordinate
    Double_t Z() const { return coord[2]; }
    /// Get at the array itself
    Double_t* GetCoord() { return coord; }

    //@}

    /// \name Overloaded operators
    /// \todo still to implement: +=, -=, *=, /=, ... ?
    //@{
    
    /// Assignment for an array.
    Coordinate& operator = (const Double_t *array)
    {
        coord[0] = array[0];
        coord[1] = array[1];
        coord[2] = array[2];
        return *this;
    }
    Coordinate& operator = (const Real_t *array)
    {
        coord[0] = array[0];
        coord[1] = array[1];
        coord[2] = array[2];
        return *this;
    }
    /// Assignment for an Coord.
    Coordinate& operator = (const Coord&q)
    {
        coord[0] = q.pos[0];
        coord[1] = q.pos[1];
        coord[2] = q.pos[2];
        return *this;
    }
    /// Asignment constructor
    Coordinate& operator = (const Coordinate q)
    {
        coord[0] = q[0];
        coord[1] = q[1];
        coord[2] = q[2];
        return *this;
    }

    /// Interface like an array.  Can set values with this, too.
    Double_t& operator [] (int i) { return coord[i]; }

    /// Const version of [] operator
    const Double_t& operator [] (int i) const { return coord[i]; }

    /// Multiply by a scalar
    Coordinate operator * (Double_t a)
    {
        return Coordinate(a * coord[0], a * coord[1], a * coord[2]);
    }

    /// Divide by a scalar
    Coordinate operator / (Double_t a)
    {
        return Coordinate(coord[0] / a, coord[1] / a, coord[2] / a);
    }

    /// Dot product -- coordinate times coordinate.
    Double_t operator * (const Coordinate& c)
    {
        return coord[0] * c.coord[0] + coord[1] * c.coord[1] +
               coord[2] * c.coord[2];
    }

    /// Adds two coordinates
    Coordinate operator + (const Coordinate& c)
    {
        return Coordinate(coord[0] + c.coord[0], coord[1] + c.coord[1],
                          coord[2] + c.coord[2]);
    }

    /// Subtracts two coordinates.
    Coordinate operator - (const Coordinate& c)
    {
        return Coordinate(coord[0] - c.coord[0], coord[1] - c.coord[1],
                          coord[2] - c.coord[2]);
    }

    void operator += (const Coordinate& c)
    {
        coord[0]+=c.coord[0];coord[1]+=c.coord[1];coord[2]+=c.coord[2];
    }
    void operator *= (const Coordinate& c)
    {
        coord[0]*=c.coord[0];coord[1]*=c.coord[1];coord[2]*=c.coord[2];
    }
    void operator *= (const Double_t& a)
    {
        coord[0]*=a;coord[1]*=a;coord[2]*=a;
    }

    /// Multiply by a scalar.
    friend Coordinate operator * (Double_t a, const Coordinate& c);

    /// Output to stream.
    friend std::ostream &operator << (std::ostream &outs, Coordinate c);
    /// Input from stream.
    friend std::istream &operator >> (std::istream &ins, Coordinate& c);

    //@}
    
    /// \name Math operations
    //@{
    /// Length of the vector (distance from origin to point).
    Double_t Length()
    {
        return sqrt(coord[0] * coord[0] + coord[1] * coord[1] + 
                    coord[2] * coord[2]);
    }

    /// Normalized vector (unit vector along direction of vector).
    Coordinate Normal()
    {
        Double_t L = Length();
        return Coordinate(coord[0] / L, coord[1] / L, coord[2] / L);
    }
    
    /// Cross product
    Coordinate Cross(Coordinate c2) {
        Coordinate c3;
        c3[0]=coord[1]*c2[2]-coord[2]*c2[1];
        c3[1]=coord[2]*c2[0]-coord[0]*c2[2];
        c3[2]=coord[0]*c2[1]-coord[1]*c2[0];
        return c3;
    }

    //@}

};

}
#endif    
