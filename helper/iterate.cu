
#include "common.h"

typedef double real_t;
using real_t = double;
using fluid_t = State<real_t, FluidStates::NSTATES,15>;

__global__ void
iterate_ADER_DG(
    const fluid_t& coeff,
    const fluid_t& state,
    const fluid_t& source,
    const fluid_t& source_grad,
    const fluid_t& flux_x,
    const fluid_t& flux_y,
    const fluid_t& flux_z,
    fluid_t& new_state)
{
  const int zoneIdx  = blockIdx.x*blockDim.x + threadIdx.x;
  const int stateIdx = blockIdx.y;
  
  constexpr int Ncoeff = fluid_t::n_coeff;
  real_t q1[Ncoeff], b[Ncoeff], x[Ncoeff];

  for (int s = 0; s < Ncoeff; s++)
    q1[s]  =  state[stateIdx](zoneIdx,s);
  for (int s = 0; s < Ncoeff; s++)
    q1[s] *= -source_grad[stateIdx](zoneIdx,s);
  for (int s = 0; s < Ncoeff; s++)
    q1[s] += source[stateIdx](zoneIdx,s);
  for (int s = 0; s < Ncoeff; s++)
    q1[s] += coeff[stateIdx](zoneIdx,s);

  for (int s = 0; s < Ncoeff; s++)
    b[s] = flux_x[stateIdx](zoneIdx,s);
#include "m2n.h" /* x = Theta_ap flux_x_pbcd */
  for (int s = 0; s < Ncoeff; s++)
    q1[s] -= x[s];
  
  for (int s = 0; s < Ncoeff; s++)
    b[s] = flux_y[stateIdx](zoneIdx,s);
#include "n2m.h" /* x = Theta_ap flux_y_apcd */
  for (int s = 0; s < Ncoeff; s++)
    q1[s] -= x[s];
  
  for (int s = 0; s < Ncoeff; s++)
    b[s] = flux_z[stateIdx](zoneIdx,(2*s+1)%Ncoeff);
#include "m2n.h" /* x = Theta_ap flux_z_abpd */
  for (int s = 0; s < Ncoeff; s++)
    q1[s] -= x[s];

  /* invert */
  for (int s = 0; s < Ncoeff; s++)
    b[s] = x[s];
#include "m2n.h"  /*  x = Labmda_ds b_abcs */
  
  for (int s = 0; s < Ncoeff; s++)
    new_state[stateIdx](zoneIdx,s) = x[s];
}

