#include <iostream>
#include <nlopt.hpp>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "common.h"

using std::get;


template<typename real_type, typename complex_type>
static auto scaled_chebyshev_basis(
    const size_t s, const size_t p,
    const real_type zmin,
    const real_type zmax,
    const std::vector<complex_type>& z)
{
  const auto m1 = 2/(zmax - zmin);
  const auto m0 = -(1 + 0.000 + zmin*m1);

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
auto optimize(const size_t p, const size_t s, const real_type beta_p, const real_type h_scale, const std::vector<complex_type> &ev_space, const bool verbose = true, const bool verbose_exception = true, const bool sanity_check = true)
{
  auto h_space = ev_space;
  real_type h_min = 0;
  real_type h_max = 0;

  for (auto& h : h_space)
  {
    using namespace std::literals;
    const auto h0 = h;
    h = h0*h_scale;
    h = h0.real()*h_scale + 1i*h0.imag()*std::pow(h_scale,beta_p);
    h_min = std::min(h_min, h.real());
  }
  assert(h_min < h_max);
  
  auto basis = scaled_chebyshev_basis(s,p,h_min,h_max,h_space);

 
  std::vector<real_type> fixed_coeff(p+1);
  real_type factorial = 1;
  for (size_t i = 0; i < p+1; i++)
  {
    factorial *= std::max(size_t{1},i);
    fixed_coeff[i] = 1.0/factorial;
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

  auto n_ineq = h_space.size();
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
  if (0)
  {
    n_ineq = 1;
    func_ineq = [](
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

      real_type fmax = -1;
      real_type res  = 0;
      int       imax = 0;
      for (size_t i = 0; i < data.npts; i++)
      {
        const auto g = 
          std::inner_product(x, x+n-1, cmat.begin() + i*n,complex_type{0});
        const auto f = real(g*conj(g))-1;
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
  }



  auto method = nlopt::LD_SLSQP;
  nlopt::opt opt(method, s+2);

  opt.set_min_objective(func,NULL);

  const real_type tol = 1.0e-13;
  opt.add_equality_mconstraint(
      func_eq, 
      &func_eq_data, 
      std::vector<real_type>(p+1,tol));
  
  std::vector<real_type> x(s+2, 1);
  opt.add_inequality_mconstraint(
      func_ineq, 
      &func_ineq_data, 
      std::vector<real_type>(n_ineq,tol)); 

  opt.set_xtol_abs(tol);
//  opt.set_ftol_abs(1.0e-14);


  real_type minf;
  decltype(opt.optimize(x,minf)) result;
  try
  {
    result = opt.optimize(x,minf);
  }
  catch (std::exception &e)
  {
    if (verbose_exception)
      printf(std::cerr, " Exception caught: % \n", e.what());
  }

  if (verbose)
  {
    printf(std::cerr, "result= % \n", result);
    printf(std::cerr, "minf= % \n", minf);
  }

  if (sanity_check)
  {
    const auto &poly = x;
    const auto &bmat = get<0>(basis);
    std::cerr << "fixed_coeff_comp= " << std::setprecision(16);
    for (int i = 0; i < p+1; i++)
    {
      const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(s+2), 0.0);
      std::cerr << res << " ";
    }
    std::cerr << std::endl;
    
    std::cerr << "fixed_coeff_gold= " ;
    for (int i = 0; i < p+1; i++)
    {
      const auto res = fixed_coeff[i];
      std::cerr << res << " ";
    }
    std::cerr << std::endl;
    
    std::cerr << "fixed_coeff_diff= " ;
    for (int i = 0; i < p+1; i++)
    {
      const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(s+2), 0.0);
      std::cerr << std::abs(res-fixed_coeff[i])/fixed_coeff[i] << " ";
    }
    std::cerr << std::endl;
    std::cerr << " ===================================== \n";


  }

  return std::make_tuple(
      std::vector<real_type>(x.begin(), x.end()-1),
      minf);
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

template<typename real_type, typename complex_type>
auto maximizeH(const size_t p, const size_t s, const real_type beta_p, const std::vector<complex_type>& ev_space)
{

  real_type h_min = 0; 
  real_type h_max = 0;
  for (auto &x : ev_space)
    h_max = std::max(h_max, std::abs(std::real(x)));
  h_max *= 2.01*s*s/p;

//  const auto max_iter = 1280;
  const auto max_steps = 1000;
  const auto tol_bisect = 0.01;

//  printf(std::cerr, "max_iter= % \n", max_iter);
  printf(std::cerr, "max_step= % \n", max_steps);
  printf(std::cerr, "tol_bisct= % \n", tol_bisect);

  bool converged = false;
  auto h = h_min;
  printf(std::cerr, " -- begin-loop -- \n");
  for (auto step = 0; step < max_steps; step++)
  {
    if ((h_max-h_min < tol_bisect*h_min) or (h_max < tol_bisect))
    {
      if (converged)
        break;
      else
        h = h_min;
    }
    else
    {
      h = 0.25*h_max + 0.75*h_min;
    }

    const auto& res = optimize(p, s, beta_p, h, ev_space, false, false,false);
    printf(std::cerr, "step= %  h_min= %  h_max= %  -- h= %  val= % \n",
        step, h_min, h_max, h, get<1>(res));

    if (std::abs(get<1>(res)) < 1.0e-12)
    {
      converged = true;
      h_min = h;
    }
    else
    {
      converged = false;
      h_max = h;
    }
  }

  assert(converged);
  const auto& res = optimize(p, s, beta_p, h, ev_space);
  printf(std::cerr, "step= %  h_min= %  h_max= %  -- h= %  val= % \n",
      -1, h_min, h_max, h, get<1>(res));

  return std::make_tuple(get<0>(res), h);
}

void test()
{
  using namespace std::literals;
  using real_type    =  double;
  using complex_type = std::complex<real_type>;

  const size_t npts = 1000;

  const real_type kappa  = 1;
  const real_type beta   = 0.5;
  const real_type beta_p = 0.25;


  const auto imag_lim = std::abs(beta);
  const auto l1 = linspace(0.0,imag_lim,npts);
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


  const auto res = std::get<0>(optimize(p, s, beta_p, h, ev_space));
  std::cout << "Coefficients: \n";
  for (auto & x : res)
  {
    std::cout << x << ", ";
  }
  std::cout << std::endl;
  {
    std::cerr << " ----------- \n";
    const auto res = std::get<0>(optimize(p, s, beta_p, h, ev_space));
    std::cout << "Coefficients: \n";
    for (auto & x : res)
    {
      std::cout << x << ", ";
    }
    std::cout << std::endl;
  }
}

void maximizeHdriver()
{
  using namespace std::literals;
  using real_type    =  double;
  using complex_type = std::complex<real_type>;

  int s = 30;
  int p = 8;

  size_t npts = 1000;
  real_type kappa  = 1;
  real_type beta   = 1; //1/10;
  real_type beta_p = 0.25;
  beta = 1.0/10;
  beta = 1.0e-14;
  beta= 0.1;

  beta_p = 0.3;
  s = 100;
//  s = std::min(s,std::max(p,static_cast<int>(std::sqrt(1.0/beta/0.15))+1));
//
  beta = 1.0;

  beta_p = 0.35;
  beta_p = 1;
  beta   = 1.0/100;

  s      = 110;
  real_type zeta = 1.0/std::sqrt(0.15);
  real_type stiff_kappa = 1000;
  s = std::max(p,static_cast<int>(zeta*std::sqrt(stiff_kappa) + 1));
#if 0 /* basic */
  beta = 1.0/stiff_kappa;
  beta_p = 1.0;
#else
  beta_p = 0.3;
  beta   = 1.0;
#endif


#if 0
  beta = 0.5;
  s    = 12;
#endif

  const auto imag_lim = std::abs(beta);
  const auto l1 = linspace(0.0,imag_lim,npts/2);
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
  
  printf(std::cerr, "npts= % \n", npts);
  printf(std::cerr, "p= % \n", p);
  printf(std::cerr, "s= % \n", s);
  printf(std::cerr, "kappa= % \n", kappa);
  printf(std::cerr, "kappa_stiff= % \n", stiff_kappa);
  printf(std::cerr, "beta= % \n",  beta);
  printf(std::cerr, "beta_p= % \n",  beta_p);


  const auto res = maximizeH<real_type,complex_type>(p,s,beta_p,ev_space);


  std::cout << "coeff = { \n";
  const auto & poly = get<0>(res);
  const auto & h    = get<1>(res);
  for (size_t i = 0; i < poly.size() -1 ; i++)
  {
    std::cout << std::setprecision(16) << poly[i] << ", ";
  }
  std::cout << poly.back()  << " };\n";
  std::cout << "h= " << h << std::endl;
  std::cout << "h/s^2= " << h/(s*s) << std::endl;
  std::cout << "h_imag= " << std::pow(h,beta_p)*beta << std::endl;
  std::cout << "h_imag/s= " << std::pow(h,beta_p)*beta/s << std::endl;
}


int main(int argc, char * argv[])
{
#if 0
  test();
#else
  maximizeHdriver();
#endif
  return 0;
}
