#define M 5
#define N 5

#include <random>
#include <chrono>
#include <iostream>
typedef double real_t;

__global__ void solve()
{
  for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M; b++)
      for (int c = 0; c <= M; c++)

  for (int idx = 0; idx < (M+1)*(M+1)*(M+1); idx += warpSize)
  {
    const int a =  idx / ((M+1)*(M+1));
    const int b =  idx /  (M+1);
    const int c = (idx %  (M+1)) + warpIdx;

    real_t y[N+1] = {0.0};
    Warp_shared<State,(M+1)*(M+1)*(M+1)> x; (qh[d*(M+1)*(M+1)*(M+1)*(M+1)]);

    for (int d = 0; d <= N; d++)
    {
      /* pre-cache qh[d][:][:][:] */
      x.load(qh[d*(M+1)*(M+1)*(M+1)*(M+1)]);

      for (int p = 0; p <= M; p++)
      {
        y[d] += matf[a*(M+1) + p]*flux_x(x[p*(M+1)*M(+1) + b*(M+1)+c]);
        y[d] += matf[b*(M+1) + p]*flux_y(x[a*(M+1)*M(+1) + p*(M+1)+c]);
        y[d] += matf[c*(M+1) + p]*flux_z(x[a*(M+1)*M(+1) + b*(M+1)+p]);
      }

      /* initial condition expansion */
      y[d] = g[a][b][c][d] - y[d];
    }

    for (int d = 0; d <= N; d++)
    {
      real_t res = 0.0;
      for (int s = 0; s <= N; s++)
      {
        res += matu[d][s]*y[s];
      }
      x[a][b][c][d] = res;
    }
  }
}

__global__ void solve3d()
{
  //  for (int a = 0; a <= M; a++)
  //   for (int b = 0; b <= M; b++)
  //    for (int c = 0; c <= M; c++)

  const int warpIdx = blockIdx.x*(blockDim.x/warpSize);
  const int laneIdx = threadIdx.x%warpSize;

  for (int idx = 0; idx < (M+1)*(M+1)*(M+1); idx < warpSize)
  {

    real_t y[N+1];
    for (int d = 0; d <= N; d++)
    {
      real_t state[3*(M+1)];
      real_t res = 0.0;
      for (int p = 0; p <= M; p++)
      {
        res += 
          matf[a][p]*flux_x(x[d][p][b][c]) + 
          matf[b][p]*flux_y(x[d][a][p][c]) + 
          matf[c][p]*flux_z(x[d][a][b][p]);
      }
      y[d] = g[a][b][c][d] - res;
    }

    for (int d = 0; d <= N; d++)
    {
      real_t res = 0.0;
      for (int s = 0; s <= N; s++)
      {
        res += matu[d][s]*y[s];
      }
      x[a][b][c][d] = res;
    }

  }
}

__global__ void solve2d()
{
  const int tx = blockDim.x*blockIdx.x + threadIdx.x;

  for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M; b++)
    {
      real_t y[N+1];
      for (int d = 0; d <= N; d++)
      {
        real_t state[3*(M+1)];
        real_t res = 0.0;
        for (int p = 0; p <= M; p++)
        {
          res += matf[a][p]*flux_x(x[d][b][p]);
          res += matf[b][p]*flux_y(x[d][p][a]) + 
        }
        y[d] = g[d][b][a] - res;
      }

      for (int d = 0; d <= N; d++)
      {
        real_t res = 0.0;
        for (int s = 0; s <= N; s++)
        {
          res += matu[d][s]*y[s];
        }
        x[d][b][a] = res;
      }
    }
}

void solve(const real_t g[M+1][M+1][M+1][N+1], real_t x[M+1][M+1][M+1][N+1], const real_t matf[M+1][M+1], const real_t matu[N+1][N+1])
{
  for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M; b++)
      for (int c = 0; c <= M; c++)
      {
        /////////////
        real_t y[N+1];
        for (int d = 0; d <= N; d++)
        {
          real_t res = 0.0;
          for (int p = 0; p <= M; p++)
          {
            res += 
              matf[a][p]*flux_x(x[p][b][c][d]) + 
              matf[b][p]*flux_y(x[a][p][c][d]) + 
              matf[c][p]*flux_z(x[a][b][p][d]);
          }
          y[d] = g[a][b][c][d] - res;
        }

        for (int d = 0; d <= N; d++)
        {
          real_t res = 0.0;
          for (int s = 0; s <= N; s++)
          {
            res += matu[d][s]*y[s];
          }
          x[a][b][c][d] = res;
        }
        //////////////
      }
}

int main(int argc, char * argv[])
{
  std::default_random_engine generator;
  std::uniform_real_distribution<real_t> distribution(0.0,1.0);

  real_t x[M+1][M+1][M+1][N+1];
  real_t g[M+1][M+1][M+1][N+1];
  real_t matf[M+1][M+1];
  real_t matu[N+1][N+1];

  for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M; b++)
      for (int c = 0; c <= M; c++)
        for (int d = 0; d <= N; d++)
        {
          x[a][b][c][d] = distribution(generator);
          g[a][b][c][d] = distribution(generator);
        }
  
  for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M; b++)
          matf[a][b] = distribution(generator);

  for (int a = 0; a <= N; a++)
    for (int b = 0; b <= N; b++)
      matu[a][b] = distribution(generator);

  int nrep= argc > 1 ? atoi(argv[1]) : 10;

  using std::cout;
  using std::endl;
  cout << "nrep= " << nrep << endl;

  cout << "Warm up ... " << endl;
  solve(g,x,matf,matu);
  solve(g,x,matf,matu);
  cout << "Benchmark ... " << endl;

  using Clock = std::chrono::high_resolution_clock;
  const auto t0 = Clock::now();
  for (int r = 0; r < nrep; r++)
    solve(g,x,matf,matu);
  const auto t1 = Clock::now();
  auto ms = std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0);
  cout << "Time per solve: " << ms.count()*1e-3/nrep<<  " us" <<  endl;


  return 0;

}
