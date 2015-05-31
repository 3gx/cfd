#include <iostream>
#include <nlopt.hpp>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "common.h"


template<typename real_type, typename complex_type>
static auto scaled_chebyshev_basis(
    const size_t s, const size_t p,
    const real_type zmin,
    const real_type zmax,
    const std::vector<complex_type>& z)
{
  const auto m1 = 2/(zmax - zmin);
  const auto m0 = -(1 + zmin*m1);

  std::vector<real_type> _b((p+1)*(s+2),0);
  {
    auto b = [&](const size_t _s, const size_t _p) -> real_type&
    {
      return _b[_p*(s+2) + _s];
    };
    b(0,0) = 1;
    b(1,0) = m0;
    b(1,1) = m1;
    for (size_t k = 0; k < s-1; k++)
    {
      b(k+2,0) = 2.0*(0 + m0*b(k+1,0)) - b(k,0);
      for (size_t j = 1; j < p+1; j++)
      {
        b(k+2,j) = 2.0*(m1*b(k+1,j-1) + m0*b(k+1,j)) - b(k,j);
      }
    }
  }

  std::vector<complex_type> _c((s+2)*z.size(),0);
  {
    auto c = [&](const size_t _s, const size_t _i) -> complex_type&
    {
      return _c[_i*(s+2) + _s];
    };
    for (size_t j = 0; j < z.size(); j++)
    {
      c(0,j) = 1;
      c(1,j) = m1*z[j] + m0;
    }
    for (size_t k = 0; k < s-1; k++)
      for (size_t i = 0; i < z.size(); i++)
      {
        c(k+2,i) = 2.0*(m1*z[i]+m0)*c(k+1,i) - c(k,i);
      }
  }

  return std::make_tuple(_b, _c);
}

template<typename real_type, typename complex_type>
auto optimize(const size_t p, const size_t s, const real_type h_scale, const std::vector<complex_type> &ev_space)
{
  auto h_space = ev_space;
  real_type h_min = 0;
  real_type h_max = 0;
  for (auto& h : h_space)
  {
    h *= h_scale;
    h_min = std::min(h_min, h.real());
  }

  assert(h_min < h_max);
  
  auto basis = scaled_chebyshev_basis(s,p,h_min,h_max,h_space);

  using std::get;
  auto &bmat = get<0>(basis);
  auto &cmat = get<1>(basis);

  struct func_eq_data_type
  {
    std::vector<real_type>* bmat_ptr;
    size_t p;
    size_t s;
  } func_eq_data;
  struct func_ineq_data_type
  {
    std::vector<complex_type> *cmat_ptr, *z_ptr;
    size_t s;
  } func_ineq_data;

  func_eq_data.bmat_ptr = &bmat;
  func_eq_data.p        = p;
  func_eq_data.s        = s;

  func_ineq_data.cmat_ptr = &cmat;
  func_ineq_data.z_ptr    = &h_space;
  func_ineq_data.s        = s;

  std::vector<real_type> fixed_coeff(p+1);
  real_type factorial = 1;
  for (size_t i = 0; i < p+1; i++)
  {
    factorial *= std::max(size_t{1},i);
    fixed_coeff[i] = 1.0/factorial;
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
  };

  nlopt::mfunc func_eq = [](
      unsigned m, double *result,
      unsigned n, const double *x,
      double *grad, /* NULL if not needed */
      void *func_data)
  {
    assert(func_data);
    const auto &data = *reinterpret_cast<func_eq_data_type*>(func_data);
    const auto& bmat = *data.bmat_ptr;
    const auto s = data.s;
    const auto p = data.p;
    assert(m == p+1);
    assert(n == s+2);
    if (grad)
    {
      std::copy(bmat.begin(), bmat.end(), grad);
    }
    for (size_t i = 0; i < m; i++)
    {
      result[i] = 
        std::inner_product(x, x+n-1, bmat.begin() + i*(s+2), 0);
    }
  };

  nlopt::mfunc func_ineq = [](
      unsigned m, double *result,
      unsigned n, const double *x,
      double *grad, /* NULL if not needed */
      void *func_data)
  {
    assert(func_data);
    const auto &data = *reinterpret_cast<func_ineq_data_type*>(func_data);
    const auto& cmat = *data.cmat_ptr;
    const auto&    z = *data.z_ptr;
    const auto s = data.s;
    assert(m == z.size());
    assert(n == s+2);
    using std::conj;
    if (grad)
    {
      for (size_t i = 0; i < m; i++)
      {
        const auto g = 
          std::inner_product(x, x+n-1, cmat.begin() + i*(s+2),complex_type{0});
        for (size_t j = 0; j < s+1; j++)
        {
          auto df = cmat[i*(s+2)+j]*conj(g);
          grad[i*(s+2)+j] = (df + conj(df)).real();
        }
        for (size_t j = s+1; j < s+2; j++)
          grad[i*(s+2)+j] = 1;

        const auto re = (g*conj(g)).real();
        result[i] = (re-1)-x[n-1];
      }
    }
    else
      for (size_t i = 0; i < m; i++)
      {
        const auto g = 
          std::inner_product(x, x+n-1, cmat.begin() + i*(s+2),complex_type{0});
        const auto re = (g*conj(g)).real();
        result[i] = (re-1)-x[n-1];
      }
  };



  nlopt::opt opt(nlopt::LD_SLSQP, s+2);

  opt.set_min_objective(func,NULL);

  const real_type tol = 1.0e-12;
  opt.add_equality_mconstraint(
      func_eq, 
      &func_eq_data, 
      std::vector<real_type>(p+1,tol));
  
  std::vector<real_type> x(s+2, 1);
  opt.add_inequality_mconstraint(
      func_ineq, 
      &func_ineq_data, 
      std::vector<real_type>(h_space.size(),tol)); 

  opt.set_xtol_rel(tol);


  real_type minf;
  const auto result = opt.optimize(x,minf);

  return x;
}

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

  size_t p = 4;
  size_t s = 30;
  real_type h = 35.58;


  const auto res = optimize(p, s, h, ev_space);
  std::cout << "Coefficients: \n";
  for (auto & x : res)
  {
    std::cout << x << ", ";
  }
  std::cout << std::endl;
}

int main(int argc, char * argv[])
{
  test();
  return 0;
}
