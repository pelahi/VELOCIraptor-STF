/*! \file Precision.h
 *  \brief this file defines several types
 */

// NBodyPrecision.h

#ifndef NBODYPRECISION_H
#define NBODYPRECISION_H

/// \typedef Double_t
/// specifies the floating point precision
/// \typedef Real_t
/// place holder for other standard precision, that is if Double_t is single precision, Real_t is double precision, and visaversa. Just for compatability issues.

/// \typedef Int_t
/// specifies the size of integers

#ifdef QUADQUADPRECISION
#include <qd/qd_real.h>
#include <qd/qd_inline.h>
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
typedef qd_real Double_t;
typedef float Real_t;
#elif QUADPRECISION
#include <qd/qd_real.h>
#include <qd/qd_inline.h>
#include <qd/dd_real.h>
#include <qd/dd_inline.h>
typedef dd_real Double_t;
typedef float Real_t;
#elif SINGLEPRECISION
typedef float Double_t;
typedef double Real_t;
#define MINVALUE 1e-16
#define MAXVALUE 1e+16
#else
#define MINVALUE 1e-32
#define MAXVALUE 1e+32
typedef double Double_t;
typedef float Real_t;
#endif

#ifdef LONGINT
typedef long long int Int_t;
typedef long long unsigned int UInt_t;
#elif SHORT
typedef short Int_t;
typedef unsigned short UInt_t;
#else
typedef int Int_t;
typedef unsigned int UInt_t;
#endif

#endif
