/*! \file endianutils.cxx
 *  \brief provides instances of global variables
 */

#include "endianutils.h"

//-- Structures and external variables

short (*BigShort) ( short s );
short (*LittleShort) ( short s );
int (*BigInt) ( int i );
int (*LittleInt) ( int i );
long long (*BigLongInt) ( long long i );
long long (*LittleLongInt) ( long long i );
float (*BigFloat) ( float f );
float (*LittleFloat) ( float f );
double (*BigDouble) ( double f );
double (*LittleDouble) ( double f );
Double_t (*BigDouble_t) ( Double_t f );
Double_t (*LittleDouble_t) ( Double_t f );

#ifdef GADGETDOUBLEPRECISION
double (*BigFLOAT) ( double f );
double (*LittleFLOAT) ( double f );
#else
float (*BigFLOAT) ( float f );
float (*LittleFLOAT) ( float f );
#endif

#ifdef GADGETSINGLEMASSPRECISION
float (*BigREAL) ( float f );
float (*LittleREAL) ( float f );
#else
double (*BigREAL) ( double f );
double (*LittleREAL) ( double f );
#endif


#ifdef RAMSESSINGLEPRECISION
float (*BigRAMSESFLOAT) ( float f );
float (*LittleRAMSESFLOAT) ( float f );
#else
double (*BigRAMSESFLOAT) ( double f );
double (*LittleRAMSESFLOAT) ( double f );
#endif

bool BigEndianSystem;
