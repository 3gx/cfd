
__global__ void matvec_glb(double xinp[], double binp[])
{
  const int N = 70;
  double x[N], b[N];

#pragma unroll
  for (int i = 0; i < N; i++)
  {
    b[i] = binp[i*blockDim.x + threadIdx.x];
  }
#ifdef N2M
#include "n2m.h"
#else
#include "m2n.h"
#endif
  
#pragma unroll
  for (int i = 0; i < N; i++)
  {
    xinp[i*blockDim.x + threadIdx.x] = x[i];
  }
}

