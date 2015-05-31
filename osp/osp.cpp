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

  std::vector<real_type> _b((s+2)*(p+1),0);
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
auto optimize(const size_t p, const size_t s, const real_type h_scale, const std::vector<complex_type> &ev_space, const bool recompute_basis = true)
{
  static auto h_space = ev_space;
  static real_type h_min = 0;
  static real_type h_max = 0;

  if (recompute_basis)
  {
    for (auto& h : h_space)
    {
      h *= h_scale;
      h_min = std::min(h_min, h.real());
    }
    assert(h_min < h_max);
  }
  
  static auto basis = scaled_chebyshev_basis(s,p,h_min,h_max,h_space);
  if (recompute_basis)
  {
    basis = scaled_chebyshev_basis(s,p,h_min,h_max,h_space);
  }

 
  static std::vector<real_type> fixed_coeff(p+1);
  if (recompute_basis)
  {
    real_type factorial = 1;
    for (size_t i = 0; i < p+1; i++)
    {
      factorial *= std::max(size_t{1},i);
      fixed_coeff[i] = 1.0/factorial;
    }
  }

  using std::get;
  struct func_eq_data_type
  {
    std::vector<real_type>* bmat_ptr, *coeff_ptr;
    size_t p;
    size_t s;
  } func_eq_data;
  func_eq_data.bmat_ptr  = &get<0>(basis);
  func_eq_data.coeff_ptr = &fixed_coeff;
  func_eq_data.p         = p;
  func_eq_data.s         = s;

  struct func_ineq_data_type
  {
    std::vector<complex_type> *cmat_ptr;
    size_t npts;
    size_t s;
  } func_ineq_data;
  func_ineq_data.cmat_ptr = &get<1>(basis);
  func_ineq_data.npts     = h_space.size();
  func_ineq_data.s        = s;

  nlopt::func func = [](
      unsigned n, const double *x,
      double *grad,
      void *func_data)
  {
    if (grad)
    {
      std::fill(grad, grad+n-1, 0);
      grad[n-1] = 1;
    }
    return x[n-1];
  };

  nlopt::mfunc func_eq = [](
      unsigned m, double *result,
      unsigned n, const double *x,
      double *grad, /* NULL if not needed */
      void *func_data)
  {
    assert(func_data);
    const auto &data  = *reinterpret_cast<func_eq_data_type*>(func_data);
    const auto& bmat  = *data.bmat_ptr;
    const auto& coeff = *data.coeff_ptr;
    const auto s = data.s;
    const auto p = data.p;
    assert(m == p+1);
    assert(n == s+2);
    assert(m*n == bmat.size());
    if (grad)
    {
      std::copy(bmat.begin(), bmat.end(), grad);
    }
    for (size_t i = 0; i < m; i++)
    {
      result[i] = -coeff[i] +
        std::inner_product(x, x+n-1, bmat.begin() + i*n, 0.0);
    }
  };

#if 1
  const auto n_ineq = h_space.size();
  nlopt::mfunc func_ineq = [](
      unsigned m, double *result,
      unsigned n, const double *x,
      double *grad, /* NULL if not needed */
      void *func_data)
  {
    assert(func_data);
    const auto &data = *reinterpret_cast<func_ineq_data_type*>(func_data);
    const auto& cmat = *data.cmat_ptr;
    const auto s = data.s;
    assert(m == data.npts);
    assert(n == s+2);
    using std::conj;
    using std::real;
    if (grad)
    {
      for (size_t i = 0; i < m; i++)
      {
        const auto g = 
          std::inner_product(x, x+n-1, cmat.begin() + i*n,complex_type{0});
        for (size_t j = 0; j < s+1; j++)
        {
          const auto df = cmat[i*n+j]*conj(g);
          grad[i*n+j] = real(df + conj(df));
        }
        grad[i*n+s+1] = -1;

        const auto re = real(g*conj(g));
        result[i] = (re-1)-x[n-1];
      }
    }
    else
    {
      for (size_t i = 0; i < m; i++)
      {
        const auto g = 
          std::inner_product(x, x+n-1, cmat.begin() + i*n,complex_type{0});
        const auto re = real(g*conj(g));
        result[i] = (re-1)-x[n-1];
      }
    }
  };
#else
  const auto n_ineq = 1;
  nlopt::mfunc func_ineq = [](
      unsigned m, double *result,
      unsigned n, const double *x,
      double *grad, /* NULL if not needed */
      void *func_data)
  {
    assert(func_data);
    const auto &data = *reinterpret_cast<func_ineq_data_type*>(func_data);
    const auto& cmat = *data.cmat_ptr;
    const auto s = data.s;
    assert(m == 1);
    assert(n == s+2);
    using std::conj;
    using std::real;

    real_type fmax = 0;
    real_type res  = 0;
    int       imax = 0;
    for (size_t i = 0; i < m; i++)
    {
      const auto g = 
        std::inner_product(x, x+n-1, cmat.begin() + i*n,complex_type{0});
      const auto f = real(g*conj(g)) - 1;
      if (f > fmax)
      {
        fmax = f;
        res =  f - x[n-1];
        imax = i;
      }
    }
    result[0] = res;

    if (grad)
    {
      const auto g = 
        std::inner_product(x, x+n-1, cmat.begin() + imax*n,complex_type{0});
      for (size_t j = 0; j < s+1; j++)
      {
        const auto df = cmat[imax*n+j]*conj(g);
        grad[j] = real(df + conj(df));
      }
      grad[s+1] = -1;
    }
  };
#endif



  auto method = nlopt::LD_SLSQP;
  nlopt::opt opt(method, s+2);

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
      std::vector<real_type>(n_ineq,tol)); 

  opt.set_xtol_rel(tol);


  real_type minf;
  decltype(opt.optimize(x,minf)) result;
  try
  {
    result = opt.optimize(x,minf);
  }
  catch (std::exception &e)
  {
    printf(std::cerr, " Exception caught: % \n", e.what());
  }

  printf(std::cerr, "result= % \n", result);
  printf(std::cerr, "minf= % \n", minf);

  return std::vector<real_type>(x.begin(), x.end()-1);
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
  const real_type beta   = 0.0;


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
      ev_space.emplace_back(-kappa + 1i*x);

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
  {
    std::cerr << " ----------- \n";
    const auto res = optimize(p, s, h, ev_space,false);
    std::cout << "Coefficients: \n";
    for (auto & x : res)
    {
      std::cout << x << ", ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char * argv[])
{
  test();
  return 0;
}
