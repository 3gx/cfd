
#include "common.h"

typedef double real_t;
using real_t = double;
using fluid_t = State<real_t, FluidStates::NSTATES,15>;


template<bool NODAL2MODAL>
__device__ void 
convert_modal2nodal2modal(const fluid_t& in, fluid_t& out)
{
  const int zoneIdx  = blockIdx.x*blockDim.x + threadIdx.x;
  const int stateIdx = blockIdx.y;

  constexpr int N = fluid_t::n_coeff;
  real_t x[N], b[N];

#pragma unroll
  for (int s = 0; s < N; s++)
  {
    b[s] = in[stateIdx](zoneIdx, s);
  }


  if (!NODAL2MODAL)
  {
#include "m2n.h"
  }
  else
  {
#include "n2m.h"
  }
  
#pragma unroll
  for (int s = 0; s < N; s++)
  {
    out[stateIdx](zoneIdx, s) = x[s];
  }
}



__global__ void 
convert_modal2nodal(const fluid_t& stateModal, fluid_t& stateNodal)
{
  convert_modal2nodal2modal<0>(stateModal, stateNodal);
}
__global__ void 
convert_nodal2modal(const fluid_t& stateNodal, fluid_t& stateModal)
{
  convert_modal2nodal2modal<1>(stateNodal, stateModal);
}
