/*! \file Random.h
 *  \brief header file for random number generators
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <cmath>
#include <Precision.h>

using namespace std;
namespace Math
{

    /*! Long period (> 2 x 10^18) random number generator of L'Ecuyer with Bays-Durham shuffle
        and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
        the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
        idum between successive deviates in a sequence. RNMX should approximate the largest floating
        value that is less than 1.
    */
    Double_t ran2(long *idum);

    ///normally distributed random number using the polar Box-Muller method that uses the \ref Math::ran2 to generate
    ///uniform random numbers
    Double_t nran2(long *idum);

}

#endif
