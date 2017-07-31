/*! \file utilities.cxx
 *  \brief this file contains an assortment of utilities 
 */

#include "TreeFrog.h"

/// Time functions
//@{
int GetMilliCount()
{
  // Something like GetTickCount but portable
  // It rolls over every ~ 12.1 days (0x100000/24/60/60)
  // Use GetMilliSpan to correct for rollover
  timeb tb;
  ftime( &tb );
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}

int GetMilliSpan( int nTimeStart )
{
  int nSpan = GetMilliCount() - nTimeStart;
  if ( nSpan < 0 )
    nSpan += 0x100000 * 1000;
  return nSpan;
}

//@}

double MyGetTime(){
#ifdef USEOPENMP
    return omp_get_wtime();
#else
    return (clock() /( (double)CLOCKS_PER_SEC));
#endif
}


