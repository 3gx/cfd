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
static void step(Param &param, SoA &f)
{
  const auto n = f.size();
  static SoA k0(n);
  static SoA k1(n);
  static SoA k2(n);
  static SoA k3(n);
  static SoA k4(n);
  static SoA k5(n);
  static SoA k6(n);

  const auto h = param.h
  auto scale = [h](SoA_value val) { return val*h; }

  for (int i = 0; i < n; i++)
  {
    k0[i] = f[i];
  }
  brusselator_rhs(param, f, k1, scale);

  for (int i = 0; i < n; i++)
  {
    k0[i]  = f[i] + Real{0.25}*k1[i];
  }
  brusselator_rhs(param, k0, k2, scale);

  for (int i = 0; i < n; i++)
  {
    k0[i]  = f[i] + Real{3.0/32.0}*k1[i] + Real{9.0/32.0}*k2[i];
  }
  brusselator_rhs(param, k0, k3, scale);

  for (int i = 0; i < n; i++)
  {
    k0[i]  = f[i] + Real{1932.0/2197.0}*k1[i] - Real{7200.0/2197.0}*k2[i] + Real{7296.0/2197}*k3[i];
  }
  brusselator_rhs(param, k0, k4, scale);

  for (int i = 0; i < n; i++)
  {
    k0[i]  = f[i] + Real{439.0/216}*k1[i] - Real{8}*k2[i] + Real{3680.0/513}*k3[i] - Real{845.0/4104}*k4[i];
  }
  brusselator_rhs(param, k0, k5, scale);

  for (int i = 0; i < n; i++)
  {
    k0[i]  = f[i] - Real{8.0/27}*k1[i] + Real{2}*k2[i] - Real{3544.0/2565}*k3[i] + Real{1859.0/4104}*k4[i] - Real{11.0/40}*k5[i];
  }
  brusselator_rhs(param, k0, k6, scale);

  auto err = Real{0};
  for (int i = 0; i < n; i++)
  {
    const auto z = f[i] + Real{16.0/135}*k1[i] + Real{6656.0/12825}*k3[i] + Real{28561.0/56430}*k4[i] - Real{9.0/50}*k5[i] + Real{2.0/55}*k6[i];
    const auto y = f[i] + Real{25.0/216}*k1[i] + Real{1408.0/ 2565}*k3[i] + Real{ 2197.0/ 4101}*k4[i] - Real{   0.2}*k5[i];
    f[i] = z;
    const auto err_loc = max(abs(z - y));
    err = std::max(err_loc, err);
  }

  const auto s = std::pow(para.eps*Real{0.5}/err, Real{0.25});
  param.h *= s;
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




