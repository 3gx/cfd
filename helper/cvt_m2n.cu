
#include "common.h"

typedef double real_t;
using real_t = double;
using fluid_t = State<real_t, FluidStates::NSTATES,15>;


__global__ void 
convert_modal2nodal(const fluid_t& stateModal, fluid_t& stateNodal)
{
  const int zoneIdx  = blockIdx.x*blockDim.x + threadIdx.x;
  const int stateIdx = blockIdx.y;

  constexpr int N = fluid_t::n_coeff;
  real_t x[N], b[N];

#pragma unroll
  for (int s = 0; s < N; s++)
  {
    b[s] = stateModal[stateIdx](zoneIdx, s);
  }

#include "m2n.h"
  
#pragma unroll
  for (int s = 0; s < N; s++)
  {
    stateNodal[stateIdx](zoneIdx, s) = x[s];
  }
}



__global__ void 
convert_nodal2modal(const fluid_t& stateNodal, fluid_t& stateModal)
{
  const int zoneIdx  = blockIdx.x*blockDim.x + threadIdx.x;
  const int stateIdx = blockIdx.y;

  constexpr int N = fluid_t::n_coeff;
  real_t x[N], b[N];

#pragma unroll
  for (int s = 0; s < N; s++)
  {
    b[s] = stateNodal[stateIdx](zoneIdx, s);
  }

#include "n2m.h"
  
#pragma unroll
  for (int s = 0; s < N; s++)
  {
    stateModal[stateIdx](zoneIdx, s) = x[s];
  }
}
