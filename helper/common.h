#pragma once

__device__ __host__ constexpr int Log2(int n, int p = 0) 
{
  return (n <= 1) ? p : Log2(n / 2, p + 1);
}
template<int N>
__device__ __host__ constexpr int Log2()
{
  static_assert((1 << Log2(N,0)) == N, "N is not a power of 2!");
  return Log2(N,0);
}

template<typename T, int Ncoeff, int VEC_SIZE>
struct Expansion
{
  static constexpr int n_coeff = Ncoeff;
  static constexpr int VEC_SIZE_LOG2 = Log2<VEC_SIZE>();

  T data[][Ncoeff][VEC_SIZE];
  __host__ __device__ T operator()(const int i /* zone */, const int s /* expansion */) const
  {
#ifdef __CUDA_ARCH__
    return __ldg(&data[i>>VEC_SIZE_LOG2][s][i & (VEC_SIZE-1)]);
#else
    return data[i>>VEC_SIZE_LOG2][s][i & (VEC_SIZE-1)];
#endif
  }
  __device__ __host__ T& operator()(const int i /* zone */, const int s /* expansion */) 
  {
    return data[i>>VEC_SIZE_LOG2][s][i & (VEC_SIZE-1)];
  }
};

template<typename T,  int Nstate, int Ncoeff, int VEC_SIZE = 32>
struct State
{
  static constexpr int n_coeff = Ncoeff;
  static constexpr int VEC_SIZE_LOG = Log2<VEC_SIZE>();

  using expansion_t = Expansion<T,Ncoeff,1<<VEC_SIZE_LOG>;
  expansion_t data[Nstate];

  __device__ __host__ const expansion_t& operator[](const int i) const { return data[i] ;}
  __device__ __host__ expansion_t& operator[](const int i) { return data[i] ;}
};

enum FluidStates
{
  DENS = 0,
  MOMX, MOMY, ETOT,
  RHOY,
  NSTATES
};
