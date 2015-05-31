#include <iostream>
#include <nlopt.hpp>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "common.h"


#if 0
template<typename real_type, typename complex_type>
std::vector<real_type>
optimize(const size_t p, const size_t s, const real_type h, const std::vector<complex_type> &ev_space)
{
  auto basis = scaled_chebyshev_basis(s,p,h_min,0,h_space);

  std::vector<real_type> fixed_coeff(p+1);
  real_type factorial = 1;
  for (size_t i = 0; i < p+1; i++)
  {
    factorial *= std::max(1,i);
    fixed_coeff[i] /= fac;
  }

  auto func = [](const std::vector<real_type> &x,
      std::vector<real_type> &grad, void*)
  {
    if (!grad.empty())
    {
      std::fill(grad.begin(), grad.end(), 0);
      grad.back() = 1;
    }
    return x.back();
  }

  auto func_eq = [&](
      std::vector<real_type> &result,
      const std::vector<real_type> &x,
      std::vector<real_type> &grad, void*)
  {
    if (!grad.empty())
    {
      std::copy(bmat.begin(), bmat.end(), grad.begin());
    }
    const auto m = result.size();
    for (size_t i = 0; i < m; i++)
    {
      result[i] = 
        std::inner_product(x.begin(), x.end()-1, bmat.begin() + i*(s+2), 0);
    }
  }

  auto func_ineq = [&](
      std::vector<real_type> &result,
      const std::vector<real_type> &x,
      std::vector<real_type> &grad, void*)
  {
    const auto m = result.size();
    if (!grad.empty())
    {
      for (size_t i = 0; i < m; i++)
      {
        const auto g = 
          std::inner_product(x.begin(), x.end()-1, cmat.begin() + i*(s+2),0);
        for (size_t j = 0; j < s+1; j++)
        {
          auto df = c[i*(s+2)+j]*std::conj(g[j]);
          grad[i*(s+2)+j] = std::real(df + std::conj(df));
        }
        for (size_t j = s+1; j < s+2; j++)
          grad[i*(s+2)+j] = 1;

        const auto re = std::real(g*std::conj(g));
        result[i] = (re-1)-x.back();
      }
    }
    else
      for (size_t i = 0; i < m; i++)
      {
        const auto g = 
          std::inner_product(x.begin(), x.end()-1, cmat.begin() + i*(s+2),0);
        const auto re = std::real(g*std::conj(g));
        result[i] = (re-1)-x.back();
      }
  }


  nlopt::opt opt(nlopt::LD_SLSQP, s+2);
  const real_type tol = 1.0e-12;
  opt.add_equality_mconstraint(func_eq, NULL, std::vector<real>(p+1,tol));
  opt.add_inequality_mconstraint(func_ineq, NULL, std::vector<real>(s+1,tol)); 

  opt.set_xtol_rel(tol);

  std::vector<real_type> x(s+2, 1);

  real_type minf;
  const auto result = opt.optmize(x,minf);

  return x;
}
#endif

template<typename real_type>
std::vector<real_type> linspace(const real_type a, const real_type b, const size_t n)
{
  std::vector<real_type> res(n);
  std::iota(res.begin(), res.end(), 0);
  for (auto& x : res)
    x = a + (b-a)*x/(n-1); 
  return res;
}

void test()
{
  using namespace std::literals;
  using real_type    =  double;
  using complex_type = std::complex<real_type>;

  const size_t npts = 1000;

  const real_type kappa  = 1;
  const real_type beta   = 0;


  const auto imag_lim = std::abs(beta);
  const auto l1 = linspace(0.0,imag_lim,50);
  const auto l2 = linspace(-kappa,0.0,npts);

  std::vector<complex_type> ev_space;
  if (imag_lim > 0)
    for (auto& x : l1)
      ev_space.emplace_back(1i*x);
  for (auto& x : l2)
    ev_space.emplace_back(1i*imag_lim + x);
  if (imag_lim > 0)
    for (auto &x : l1)
      ev_space.emplace_back(1i*x - kappa);


  printf(std::cerr, "ev_space.size()= % \n" , ev_space.size());
}

int main(int argc, char * argv[])
{
  test();
  return 0;
}
