#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "printf.h"
#include "zip_iterator.h"

template<typename Real>
struct Params
{
  using value_type = Real;
  Real time;
  Real cfl;
  Real dx;
  Real alpha;

  real_t dt() const {return cfl * diff/square(dx);}
};

template<typename Param, typename SoA>
static void brusselator_rhs(const Param &param, SoA &f, SoA &rhs)
{
  using std::get;
  using Real = typename SoA::value_type;
  const auto n = f.size();

  /* BC */
  auto u = [](const int i) { return get<0>(f[i]); }
  auto v = [](const int i) { return get<1>(f[i]); }
  u[0] = u[n-1] = 1;
  v[0] = v[n-1] = 3;

  
  const auto c= param.alpha/square(param.dx);
  for (int i = 1; i < n-1; i++)
  {
    rhs[i] = {
      1 - 4*u[i] + square(u[i])*v[i] + c * (u[i-1] - 2*u[i] + u[i+1]),
      0 + 3*u[i] - square(u[i])*v[i] + c * (v[i-1] - 2*v[i] + v[i+1])
    }
  }
  rhs[0] = rhs[n-1] = {Real{0},Real{0}};
}

template<typename Param, typename SoA>
static void brusselator_ic(const Param &param, SoA &f)
{
  using std::get;
  using Real = typename SoA::value_type;
  const auto n = f.size();
  const auto pi = std::acos(Real{-1});

  for (int i = 0; i < n; i++)
  {
    f[i] = {
      1 + std::sin(2*pi*i/(n-1));
      Real{3}
    };
  }

}

template<typename Param, typename SoA>
void dump2file(const Param &param, const SoA &f, const std::string fileName)
{
  using std::get;
  using Real = typename SoA::value_type;
  const int n = f.size();
  std::ofstream fout(fileName);
  fout << "# time= " << std::setprecision(16) << param.time << std::endl;
  for (int i = 0; i < n; i++)
  {
    const auto x = Real{i}/(n-1);
    fout << set::precision(16) << x << "\t" << get<0>(f[i]) << "\t" <<  get<1>(f[i]) << std::endl;
  }
}

int main(int argc, char * argv[])
{
  const int ncell = argc > 1 ? atoi(argv[1]) : 100;
  fprintf(std::cerr, "ncell= %\n", ncell);

  const int niter = argc > 2 ? atoi(argv[2]) : 10;
  fprintf(std::cerr, "niter= %\n", niter);
  

  using real_t = double;

  using params_t = Params<real_t>;
  using vector_t = std::vector<real_t>;

  auto params = params_t{};
  params.dx   = 1;
  params.diff = 1;
  params.cfl  = 0.50;  /* stable for cfl <= 0.5 */
//  params.cfl  = 0.252;

  vector_t f(ncell+2);

  set_ic(params,f);
  dump2file(f,"ic.txt");

  for (int iter = 0; iter < niter; iter++)
  {
    fprintf(std::cerr, "iter= %\n", iter);
#if 1
    compute_update<FreeBC>(params,f);
//    compute_update<PeriodicBC>(params,f);
#elif 1
    compute_update_dg1<FreeBC>(FreeBC,params,f);
#endif
    dump2file(f,"iter"+std::to_string(iter));
  }




  return 0;  

}




