///Math code
#include <NBodyMath.h>
using namespace Math;

#ifndef ENDIANUTILS_H
#define ENDIANUTILS_H

//functions that reverse endian, check endian and some useful pointers. This code is from 
//http://www.gamedev.net/reference/articles/article2091.asp
//I've also copied some of their comments. 
//255 in hex is 0x00FF
inline short ShortSwap( short s )
{
/*  unsigned char b1, b2;
  b1 = s & 255;
  b2 = ( s>>8 ) & 255;
  return (b1 << 8) + b2;*/
  int size=2;
  union
  {
     short value;
     unsigned char b[2];
  } dat1, dat2;
  dat1.value = s;
  for (int j = 0; j < size / 2; ++j)
  {
     dat2.b[j] = dat1.b[size - 1 - j];
     dat2.b[size - 1 - j] = dat1.b[j];
  }
  return dat2.value;
}

//This function is fairly straightforward once you wrap your head around the bit math. We take apart the two bytes of the short parameter s with some simple bit math and then glue them back together in reverse order. If you understand bit shifts and bit ANDs, this should make perfect sense. As a companion to ShortSwap, we'll have ShortNoSwap, which is very simple:

inline short ShortNoSwap( short s ){return s;}

//Next, we want to swap 4 byte integers:
inline int IntSwap (int i)
{
/*  unsigned char b1, b2, b3, b4;
  b1 = i & 255;
  b2 = ( i>>8 ) & 255;
  b3 = ( i>>16 ) & 255;
  b4 = ( i>>24 ) & 255;
  return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;*/
  int size=4;
  union
  {
     int value;
     unsigned char b[4];
  } dat1, dat2;
  dat1.value = i;
  for (int j = 0; j < size / 2; ++j)
  {
     dat2.b[j] = dat1.b[size - 1 - j];
     dat2.b[size - 1 - j] = dat1.b[j];
  }
  return dat2.value;
}
inline int IntNoSwap( int i ){return i;}

//Next, we want to swap 8 byte integers:
inline long long LongIntSwap (long long i)
{
  int size=8;
  union
  {
     long long int value;
     unsigned char b[8];
  } dat1, dat2;
  dat1.value = i;
  for (int j = 0; j < size / 2; ++j)
  {
     dat2.b[j] = dat1.b[size - 1 - j];
     dat2.b[size - 1 - j] = dat1.b[j];
  }
  return dat2.value;
}
inline long long LongIntNoSwap( long long i ){return i;}

//LongSwap is more or less the same idea as ShortSwap, but it switches around 4 bytes instead of 2. Again, this is straightforward bit math.
//Lastly, we need to be able to swap floats:

inline float FloatSwap( float f )
{
  int size=4;
  union
  {
    float f;
    unsigned char b[4];
  } dat1, dat2;

  dat1.f = f;
  for (int j = 0; j < size / 2; ++j)
  {
     dat2.b[j] = dat1.b[size - 1 - j];
     dat2.b[size - 1 - j] = dat1.b[j];
  }
/*  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];*/
  return dat2.f;
}
inline float FloatNoSwap( float f ){return f;}
inline double DoubleSwap( double f )
{
  int size=8;
  union
  {
    double f;
    unsigned char b[8];
  } dat1, dat2;

  dat1.f = f;
  for (int j = 0; j < size / 2; ++j)
  {
     dat2.b[j] = dat1.b[size - 1 - j];
     dat2.b[size - 1 - j] = dat1.b[j];
  }
/*  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];*/
  return dat2.f;
}
inline double DoubleNoSwap( double f ){return f;}

//This is code I've written which I think will do a generic switch. Use a template
/*template <class T, unsigned int size> inline T SwapValue(T value)
{
  union
  {
     T value;
     char bytes[size];
  } in, out;

  in.value = value;

  for (unsigned int i = 0; i < size / 2; ++i)
  {
     out.bytes[i] = in.bytes[size - 1 - i];
     out.bytes[size - 1 - i] = in.bytes[i];
  }

  return out.value;
}
template <class T, unsigned int size> inline T SwapValue(T value){return vale;}*/

//The next part of our implementation is the clever twist
//that defines this method of endian independence.
//We're going to use function pointers to automatically select the correct endian for us.
//We'll put their declarations in a header file

extern short (*BigShort) ( short s );
extern short (*LittleShort) ( short s );
extern int (*BigInt) ( int i );
extern int (*LittleInt) ( int i );
extern long long (*BigLongInt) ( long long i );
extern long long (*LittleLongInt) ( long long i );
extern float (*BigFloat) ( float f );
extern float (*LittleFloat) ( float f );
extern double (*BigDouble) ( double f );
extern double (*LittleDouble) ( double f );
extern Double_t (*BigDouble_t) ( Double_t f );
extern Double_t (*LittleDouble_t) ( Double_t f );

#ifdef GADGETDOUBLEPRECISION
extern double (*BigFLOAT) ( double f );
extern double (*LittleFLOAT) ( double f );
#else
extern float (*BigFLOAT) ( float f );
extern float (*LittleFLOAT) ( float f );
#endif

#ifdef GADGETSINGLEMASSPRECISION
extern float (*BigREAL) ( float f );
extern float (*LittleREAL) ( float f );
#else
extern double (*BigREAL) ( double f );
extern double (*LittleREAL) ( double f );
#endif

#ifdef RAMSESSINGLEPRECISION
extern float (*BigRAMSESFLOAT) ( float f );
extern float (*LittleRAMSESFLOAT) ( float f );
#else
extern double (*BigRAMSESFLOAT) ( double f );
extern double (*LittleRAMSESFLOAT) ( double f );
#endif


//In order to initialize all of these function pointers,
//we'll need a function to detect the endian and set them.
//This function will make use of the byte array cast
//similar to earlier example of a byte cast that assumes a certain endian.

extern bool BigEndianSystem;  //you might want to extern this

inline void InitEndian( void )
{
  unsigned char SwapTest[2] = { 1, 0 };
  //bool BigEndianSystem;

  //if( *(short *) SwapTest == 0 )
  if( *(short *) SwapTest == 1 )
  {
    //little endian
    BigEndianSystem = false;

    //set func pointers to correct funcs
    BigShort = ShortSwap;
    LittleShort = ShortNoSwap;
    BigInt = IntSwap;
    LittleInt = IntNoSwap;
    BigLongInt = LongIntSwap;
    LittleLongInt = LongIntNoSwap;
    BigFloat = FloatSwap;
    LittleFloat = FloatNoSwap;
    BigDouble = DoubleSwap;
    LittleDouble = DoubleNoSwap;
#ifdef SINGLEPRECISION
    LittleDouble_t=FloatNoSwap;
    BigDouble_t=FloatSwap;
#else
    LittleDouble_t=DoubleNoSwap;
    BigDouble_t=DoubleSwap;
#endif

#ifdef GADGETDOUBLEPRECISION
    BigFLOAT= DoubleSwap;
    LittleFLOAT = DoubleNoSwap;
#else
    BigFLOAT= FloatSwap;
    LittleFLOAT = FloatNoSwap;
#endif
#ifdef GADGETSINGLEMASSPRECISION
    BigREAL= FloatSwap;
    LittleREAL = FloatNoSwap;
#else
    BigREAL= DoubleSwap;
    LittleREAL = DoubleNoSwap;
#endif

#ifdef RAMSESSINGLEPRECISION
    BigRAMSESFLOAT= FloatSwap;
    LittleRAMSESFLOAT = FloatNoSwap;
#else
    BigRAMSESFLOAT= DoubleSwap;
    LittleRAMSESFLOAT = DoubleNoSwap;
#endif
      
    }
  else
  {
    //big endian
    BigEndianSystem = true;

    BigShort = ShortNoSwap;
    LittleShort = ShortSwap;
    BigInt = IntNoSwap;
    LittleInt = IntSwap;
    BigLongInt = LongIntNoSwap;
    LittleLongInt = LongIntSwap;
    BigFloat = FloatNoSwap;
    LittleFloat = FloatSwap;
    BigDouble = DoubleNoSwap;
    LittleDouble = DoubleSwap;
#ifdef SINGLEPRECISION
    LittleDouble_t=FloatSwap;
    BigDouble_t=FloatNoSwap;
#else
    LittleDouble_t=DoubleSwap;
    BigDouble_t=DoubleNoSwap;
#endif
#ifdef GADGETDOUBLEPRECISION
    BigFLOAT= DoubleNoSwap;
    LittleFLOAT = DoubleSwap;
#else
    BigFLOAT= FloatNoSwap;
    LittleFLOAT = FloatSwap;
#endif
#ifdef GADGETSINGLEMASSPRECISION
    BigREAL= FloatNoSwap;
    LittleREAL = FloatSwap;
#else
    BigREAL= DoubleNoSwap;
    LittleREAL = DoubleSwap;
#endif
#ifdef RAMSESSINGLEPRECISION
    BigRAMSESFLOAT= FloatNoSwap;
    LittleRAMSESFLOAT = FloatSwap;
#else
    BigRAMSESFLOAT= DoubleNoSwap;
    LittleRAMSESFLOAT = DoubleSwap;
#endif
  }
}

//now with this code, I can alter how structures are read and make it endian independent.

#endif
