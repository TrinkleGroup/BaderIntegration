#ifndef __DCOMP_H
#define __DCOMP_H
#include <math.h>

/*
  Program: dcomp.H
  Author:  D. Trinkle
  Date:    Sept. 27, 2002
  Purpose: Define dcomp, zero, and round (obsolete in Mac OS X -- math.h)
*/

// Some general functions, etc.

const double TOLER = 1e-7;

inline int dcomp (double a, double b) 
{
  if (a>b)
    return ( (a-b) < TOLER );
  else
    return ( (b-a) < TOLER );
}

inline int zero (double a) 
{
  return ( (a < TOLER) && (a > -TOLER) );
}

/*
inline int round (double x) 
{
  if (x > 0)
    return (int) (x+0.5);
  else
    return -((int) (-x+0.5));
}
*/

#endif
