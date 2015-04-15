
#include "common.h"

typedef double real_t;
using real_t = double;
using fluid_t = State<real_t, FluidStates::NSTATES,15>;

#include <utility>
template< bool B, class T = void >
using enableIf = typename std::enable_if<B,T>::type;

template<int K, int J, int N> __forceinline__ __device__ 
real_t decompLUimpl3(const real_t D[N*N])
{
  real_t sum = 0;
#pragma unroll
  for (int p = 0; p < K; p++)
    sum += D[K*N+p]*D[p*N+J];
  return sum;
}
template<int K, int N, int J> __forceinline__ __device__ 
auto decompLUimpl1(const real_t S[N*N], real_t D[N*N]) -> enableIf<(J==N)> {}
template<int K, int N, int J = K> __forceinline__ __device__ 
auto decompLUimpl1(const real_t S[N*N], real_t D[N*N]) -> enableIf<(J<N)>
{
  const auto sum = decompLUimpl3<K,J,N>(D);
  D[K*N+J] = S[K*N+J]-sum;
  decompLUimpl1<K,N,J+1>(S,D);
}

template<int K, int N, int I> __forceinline__ __device__ 
auto decompLUimpl2(const real_t S[N*N], real_t D[N*N]) -> enableIf<(I==N)> {}
template<int K, int N, int I> __forceinline__ __device__ 
auto decompLUimpl2(const real_t S[N*N], real_t D[N*N]) -> enableIf<(I<N)>
{
  const auto sum = decompLUimpl3<I,K,N>(D);
  D[I*N+K] = (S[I*N+K] - sum)*(real_t(1.0)/D[K*N+K]);
  decompLUimpl2<K,N,I+1>(S,D);
}

template<int N, int K> __forceinline__ __device__ 
auto decompLUimpl(const real_t S[N*N], real_t D[N*N]) -> enableIf<K==N> {}
template<int N, int K = 0> __forceinline__ __device__ 
auto decompLUimpl(const real_t S[N*N], real_t D[N*N]) -> enableIf<(K<N)>
{
  decompLUimpl1<K,N,K>(S,D);
  decompLUimpl2<K,N,K+1>(S,D);
  decompLUimpl<N,K+1>(S,D);
}

template<int N> __forceinline__ __device__
void decompLU(const real_t S[N*N], real_t D[N*N])  
{
  decompLUimpl<N>(S,D);
}

template<int N>
__forceinline__ __device__
void decompLU1(const real_t S[N*N], real_t D[N*N])  
{
#pragma unroll
  for(int k = 0; k < N; k++)
  {
    // lvl1
#pragma unroll
    for(int j = k; j < N; j++)
    {
      // lvl3
      real_t sum = 0;
#pragma unroll
      for(int p = 0; p < k; p++)
        sum += D[k*N+p]*D[p*N+j];
      D[k*N+j] = (S[k*N+j]-sum); // not dividing by diagonals
    }
     
    // lvl2
#pragma unroll
    for(int i = k+1; i < N; i++)
    {
      real_t sum=0.;
#pragma unroll
      for(int p = 0; p < k; p++)
        sum += D[i*N+p]*D[p*N+k];
      D[i*N+k] = (S[i*N+k]-sum)*(real_t(1.0)/D[k*N+k]);
    }
  }
}

template<int N, int I = 0> __forceinline__ __device__
auto solveLUimpl1(const real_t LU[N*N], real_t b[N]) -> enableIf<(I==N)> {}
template<int N, int I = 0> __forceinline__ __device__
auto solveLUimpl1(const real_t LU[N*N], real_t b[N]) -> enableIf<(I<N)>
{
  real_t sum = 0;
#pragma unroll
  for (int k = 0; k < I; k++)
    sum = LU[I*N+k]*b[k];
  b[I] = b[I] - sum;
  solveLUimpl1<N,I+1>(LU,b);
}

template<int N, int I = 0> __forceinline__ __device__
auto solveLUimpl2(const real_t LU[N*N], real_t x[N]) -> enableIf<(I==-1)> {}
template<int N, int I = N-1> __forceinline__ __device__
auto solveLUimpl2(const real_t LU[N*N], real_t x[N]) -> enableIf<(I>=0)>
{
  real_t sum = 0;
#pragma unroll
  for (int k = 0; k < I; k++)
    sum = LU[I*N+k]*x[k];
  x[I] = (x[I] - sum)*(real_t(1.0)/LU[I*N+I]);
  solveLUimpl2<N,I-1>(LU,x);
}
template<int N> __forceinline__ __device__
void solveLU(const real_t LU[N*N], const real_t b[N], real_t x[N])
{
#pragma unroll
  for (int i = 0; i < N; i++)
    x[i] = b[i];
  solveLUimpl1<N>(LU,x);
  solveLUimpl2<N>(LU,x);
}

template<int N>
__forceinline__ __device__
void solveLU1(const real_t LU[N*N], const real_t b[N], real_t x[N])
{
   real_t y[N];
#pragma unroll
   for (int i = 0; i < N; i++)
   {
      real_t sum = 0;
#pragma unroll
      for (int k = 0 ;k < i; k++)
        sum += LU[i*N+k]*y[k];
      y[i] = (b[i]-sum); // not dividing by diagonals
   }

#pragma unroll
   for (int i = N-1; i >= 0; i--)
   {
      real_t sum = 0;
#pragma unroll
      for (int k = i+1; k < N; k++)
        sum += LU[i*N+k]*x[k];
      x[i]=(y[i]-sum)*(real_t(1.0)/LU[i*N+i]);
   }
}

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
  constexpr bool SOURCE=true;
  const int zoneIdx  = blockIdx.x*blockDim.x + threadIdx.x;
  const int stateIdx = blockIdx.y;
  
  constexpr int Ncoeff = fluid_t::n_coeff;
  constexpr int Morder = 3;
  real_t q1[Ncoeff], b[Ncoeff], x[Ncoeff];

  if (SOURCE)
  {
    for (int s = 0; s < Ncoeff; s++)
      q1[s]  =  state[stateIdx](zoneIdx,s);
    for (int s = 0; s < Ncoeff; s++)
      q1[s] *= -source_grad[stateIdx](zoneIdx,s);
    for (int s = 0; s < Ncoeff; s++)
      q1[s] += source[stateIdx](zoneIdx,s);
  }
  else
  {
    for (int s = 0; s < Ncoeff; s++)
      q1[s] = real_t(0.0);
  }
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
    b[s] = x[s];

  real_t Lambda[Morder*Morder];
  real_t LU    [Morder*Morder];

#pragma unroll
  for (int s = 0; s < Ncoeff*Ncoeff; s++)
    if (s < Morder*Morder)
      Lambda[s] = b[s%Ncoeff] + x[s%Ncoeff];

  decompLU<Morder>(Lambda,LU);
  solveLU<Morder>(LU,b,x);




  
  for (int s = 0; s < Ncoeff; s++)
    new_state[stateIdx](zoneIdx,s) = x[s];
}

