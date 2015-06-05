#include <sys/time.h>
#pragma once

inline double rtc(void)
{
  struct timeval Tvalue;
  double etime;
  struct timezone dummy;

  //gettimeofday(&Tvalue,NULL);
  gettimeofday(&Tvalue,&dummy);
  etime =  (double) Tvalue.tv_sec +
    1.e-6*((double) Tvalue.tv_usec);
  return etime;
}

inline void rdtsc_(unsigned long long * count)
{
  unsigned int low,high;
  __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
  *count =  low + (((long) high)<<32);
}

double rtc_(void)
{
  struct timeval Tvalue;
  double etime;
  struct timezone dummy;

  //gettimeofday(&Tvalue,NULL);
  gettimeofday(&Tvalue,&dummy);
  etime =    (double) Tvalue.tv_sec
    + 1.e-6*((double) Tvalue.tv_usec);
  return etime;
}
