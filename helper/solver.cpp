#define M 5
#define N 5

#include <random>
#include <chrono>
#include <iostream>
typedef double real_t;

void solve(const real_t g[M+1][M+1][M+1][N+1], real_t x[M+1][M+1][M+1][N+1], const real_t matf[M+1][M+1], const real_t matu[N+1][N+1])
{
  for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M; b++)
      for (int c = 0; c <= M; c++)
      {
        real_t y[N+1];
        for (int d = 0; d <= N; d++)
        {
          real_t res = 0.0;
          for (int p = 0; p <= M; p++)
          {
            res += 
              matf[a][p]*x[p][b][c][d] + 
              matf[b][p]*x[a][p][c][d] + 
              matf[c][p]*x[a][b][p][d];
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
  cout << "Time per solve: " << ms.count()*1.0/nrep<< endl;


  return 0;

}
