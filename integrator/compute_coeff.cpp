#include <iostream>
#include <nlopt.hpp>
#include <complex>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cassert>
#include "common.h"

template<typename Real, typename Complex>
class OptimizerT
{
  using real_type = Real;
  using complex_type = Complex;
  private:
    size_t _s;   /* number of stages */
    size_t _p;   /* order            */
    std::vector<real_type> _bmat;
    std::vector<complex_type> _cmat;
    std::vector<complex_type> _z, _ev_space;
    real_type _zmin;

    std::vector<real_type> _solution;
    real_type _fmin;

    real_type _tol;
    bool _verbose;

  public:
    OptimizerT(size_t s, size_t p, std::vector<complex_type> ev_space) :
      _s(s), _p(p), 
      _ev_space(std::move(ev_space)),
      _tol(1.0e-14),
      _verbose(true)
    {
      rescale_ev(1,1);
    }
    void set_s(const size_t s) { _s = s; }
    auto get_s() const { return _s; }

    void set_p(const size_t p) { _p = p; }
    auto get_p() const { return _p; }

    void set_tol(const real_type tol) { _tol = tol; }
    auto get_tol() const { return _tol; }

    void set_verbose() { _verbose = true; }
    void unset_verbose() { _verbose = false; }

    auto get_zmin() const { return _zmin; }

    auto get_solution() const { return std::vector<real_type>(_solution.begin(), _solution.end()-1); }
    auto get_fmin() const { return _fmin; }


  private:

    void rescale_ev(const real_type h_real, const real_type h_imag)
    {
      using namespace std::literals;
      using std::real;
      using std::imag;
      _z.clear();
      _zmin = 0;
      for (auto ev : _ev_space)
      {
        const auto z = real(ev)*h_real + 1i*imag(ev)*h_imag;
        _zmin = std::min(_zmin, real(z));
        _z.push_back(z);
      }
      assert(_zmin < 0);
      assert(_z.size() == _ev_space.size());
    }

    void scaled_chebyshev_basis()
    {
      const real_type zmax = 0;
      const auto zmin = _zmin;
      const auto& z = _z;

      const auto m1 = 2/(zmax - zmin);
      const auto m0 = -(1 + zmin*m1);

      _bmat.resize((_s+2)*(_p+1));
      std::fill(_bmat.begin(), _bmat.end(), 0);

      {
        auto b = [&](const size_t s, const size_t p) -> real_type&
        {
          return _bmat[p*(_s+2) + s];
        };
        b(0,0) = 1;
        b(1,0) = m0;
        b(1,1) = m1;
        for (size_t k = 0; k < _s-1; k++)
        {
          b(k+2,0) = 2.0*(0 + m0*b(k+1,0)) - b(k,0);
          for (size_t j = 1; j < _p+1; j++)
          {
            b(k+2,j) = 2.0*(m1*b(k+1,j-1) + m0*b(k+1,j)) - b(k,j);
          }
        }
      }

      _cmat.resize((_s+2)*z.size());
      std::fill(_cmat.begin(), _cmat.end(), 0);
      {
        auto c = [&](const size_t s, const size_t i) -> complex_type&
        {
          return _cmat[i*(_s+2) + s];
        };
        for (size_t j = 0; j < z.size(); j++)
        {
          c(0,j) = 1;
          c(1,j) = m1*z[j] + m0;
        }
        for (size_t k = 0; k < _s-1; k++)
          for (size_t i = 0; i < z.size(); i++)
          {
            c(k+2,i) = 2.0*(m1*z[i]+m0)*c(k+1,i) - c(k,i);
          }
      }
    }

  public:

    auto optimize(const real_type h_real , const real_type h_imag = 1, const std::vector<real_type> x_guess = std::vector<real_type>())
    {
      rescale_ev(h_real, h_imag);
      scaled_chebyshev_basis();

      std::vector<real_type> fixed_coeff(_p+1);
      real_type factorial = 1;
      for (size_t i = 0; i < _p+1; i++)
      {
        factorial *= std::max(size_t{1},i);
        fixed_coeff[i] = 1.0/factorial;
      }

      struct func_eq_data_type
      {
        std::vector<real_type>* bmat_ptr, *coeff_ptr;
        size_t p;
        size_t s;
      } func_eq_data;
      func_eq_data.bmat_ptr  = &_bmat;
      func_eq_data.coeff_ptr = &fixed_coeff;
      func_eq_data.p         = _p;
      func_eq_data.s         = _s;

      struct func_ineq_data_type
      {
        std::vector<complex_type> *cmat_ptr;
        size_t npts;
        size_t s;
      } func_ineq_data;
      func_ineq_data.cmat_ptr = &_cmat;
      func_ineq_data.npts     = _z.size();
      func_ineq_data.s        = _s;

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
          result[i] = -coeff[i] + std::inner_product(x, x+n-1, bmat.begin() + i*n, 0.0);
        }
      };

      auto n_ineq = _z.size();
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
  
      auto method = nlopt::LD_SLSQP;
      nlopt::opt opt(method, _s+2);

      opt.set_min_objective(func,NULL);

      const real_type tol = _tol;
      opt.add_equality_mconstraint(
          func_eq, 
          &func_eq_data, 
          std::vector<real_type>(_p+1,tol));

      std::vector<real_type> &x = _solution;
      x.resize(_s+2, 1);

      if (x_guess.empty())
      {
        std::generate(x.begin(), x.end(), [](){ return drand48(); });
//        std::fill(x.begin(), x.end(), 1);
      }
      else
      {
        std::copy(x_guess.begin(), x_guess.end(), x.begin());
        x.back() = 0;
      }
      opt.add_inequality_mconstraint(
          func_ineq, 
          &func_ineq_data, 
          std::vector<real_type>(n_ineq,tol)); 

//      opt.set_xtol_rel(1.0e-20);
      opt.set_xtol_abs(tol);
//     opt.set_ftol_abs(tol);
//      opt.set_ftol_rel(tol);

      real_type& minf = _fmin;
      decltype(opt.optimize(x,minf)) result;
      try
      {
        result = opt.optimize(x,minf);
      }
      catch (std::exception &e)
      {
        if (_verbose)
          printf(std::cerr, " Exception caught: % \n", e.what());
      }

      if (_verbose)
      {
//        printf(std::cerr, "result= % \n", result);
        printf(std::cerr, "minf= % \n", minf);
      }

      if (_verbose)
      {
        const auto &poly = x;
        const auto &bmat = _bmat;
#if 0
        std::cerr << "fixed_coeff_comp= " << std::setprecision(16);
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          std::cerr << res << " ";
        }
        std::cerr << std::endl;

        std::cerr << "fixed_coeff_gold= " ;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = fixed_coeff[i];
          std::cerr << res << " ";
        }
        std::cerr << std::endl;
#endif
        std::cerr << "fixed_coeff_diff= " ;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          std::cerr << std::abs(res-fixed_coeff[i])/fixed_coeff[i] << " ";
        }
        std::cerr << std::endl;
//        std::cerr << " ===================================== \n";
      }
      real_type norm = 0;
      {
        const auto &poly = x;
        const auto &bmat = _bmat;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          const auto f = (res-fixed_coeff[i])/fixed_coeff[i];
          norm += f*f;
        }
        norm = std::sqrt(norm/(_p+1));
      }
      return norm;

    }
};

template<typename Real, typename Complex>
class OptimizerL_CMP
{
  using real_type = Real;
  using complex_type = Complex;
  private:
    size_t _s;   /* number of stages */
    size_t _p;   /* order            */
    std::vector<real_type> _bmat;
    std::vector<complex_type> _cmat;
    std::vector<complex_type> _z, _ev_space;
    real_type _zmin;

    std::vector<real_type> _solution;
    real_type _fmin;

    real_type _tol;
    bool _verbose;

  public:
    OptimizerL_CMP(size_t s, size_t p, std::vector<complex_type> ev_space) :
      _s(s), _p(p), 
      _ev_space(std::move(ev_space)),
      _tol(1.0e-13),
      _verbose(true)
    {
      rescale_ev(1,1);
    }
    void set_s(const size_t s) { _s = s; }
    auto get_s() const { return _s; }

    void set_p(const size_t p) { _p = p; }
    auto get_p() const { return _p; }

    void set_tol(const real_type tol) { _tol = tol; }
    auto get_tol() const { return _tol; }

    void set_verbose() { _verbose = true; }
    void unset_verbose() { _verbose = false; }

    auto get_zmin() const { return _zmin; }

    auto get_solution() const { return std::vector<real_type>(_solution.begin(), _solution.end()-1); }
    auto get_fmin() const { return _fmin; }


  private:

    void rescale_ev(const real_type h_real, const real_type h_imag)
    {
      using namespace std::literals;
      using std::real;
      using std::imag;
      _z.clear();
      _zmin = 0;
      for (auto ev : _ev_space)
      {
        const auto z = real(ev)*h_real + 1i*imag(ev)*h_imag;
        _zmin = std::min(_zmin, real(z));
        _z.push_back(z);
      }
      assert(_zmin < 0);
      assert(_z.size() == _ev_space.size());
    }

    void scaled_legendre_basis()
    {
      const real_type zmax = 0;
      const auto zmin = _zmin;
      const auto& z = _z;

      const auto m1 = 2/(zmax - zmin);
      const auto m0 = -(1 + zmin*m1);

      _bmat.resize((_s+2)*(_p+1));
      std::fill(_bmat.begin(), _bmat.end(), 0);

      {
        auto b = [&](const size_t s, const size_t p) -> real_type&
        {
          return _bmat[p*(_s+2) + s];
        };
        b(0,0) = 1;
        b(1,0) = m0;
        b(1,1) = m1;
        for (size_t k = 0; k < _s-1; k++)
        {
          const auto a_coeff = ( 2.0*(k+2)-1) /(k+2);
          const auto b_coeff = (-1.0*(k+2 -1))/(k+2);
          b(k+2,0) = a_coeff*(0 + m0*b(k+1,0)) + b_coeff*b(k,0);
          for (size_t j = 1; j < _p+1; j++)
          {
            b(k+2,j) = a_coeff*(m1*b(k+1,j-1) + m0*b(k+1,j)) + b_coeff*b(k,j);
          }
        }
      }

      _cmat.resize((_s+2)*z.size());
      std::fill(_cmat.begin(), _cmat.end(), 0);
      {
        auto c = [&](const size_t s, const size_t i) -> complex_type&
        {
          return _cmat[i*(_s+2) + s];
        };
        for (size_t j = 0; j < z.size(); j++)
        {
          c(0,j) = 1;
          c(1,j) = m1*z[j] + m0;
        }
        for (size_t k = 0; k < _s-1; k++)
        {
          const auto a_coeff = ( 2.0*(k+2)-1) /(k+2);
          const auto b_coeff = (-1.0*(k+2 -1))/(k+2);
          for (size_t i = 0; i < z.size(); i++)
          {
            c(k+2,i) = a_coeff*(m1*z[i]+m0)*c(k+1,i) + b_coeff*c(k,i);
          }
        }
      }
    }

  public:

    auto optimize(const real_type h_real , const real_type h_imag = 1, const std::vector<real_type> x_guess = std::vector<real_type>())
    {
      rescale_ev(h_real, h_imag);
      scaled_legendre_basis();

      std::vector<real_type> fixed_coeff(_p+1);
      real_type factorial = 1;
      for (size_t i = 0; i < _p+1; i++)
      {
        factorial *= std::max(size_t{1},i);
        fixed_coeff[i] = 1.0/factorial;
      }

      struct func_eq_data_type
      {
        std::vector<real_type>* bmat_ptr, *coeff_ptr;
        size_t p;
        size_t s;
      } func_eq_data;
      func_eq_data.bmat_ptr  = &_bmat;
      func_eq_data.coeff_ptr = &fixed_coeff;
      func_eq_data.p         = _p;
      func_eq_data.s         = _s;

      struct func_ineq_data_type
      {
        std::vector<complex_type> *cmat_ptr;
        size_t npts;
        size_t s;
      } func_ineq_data;
      func_ineq_data.cmat_ptr = &_cmat;
      func_ineq_data.npts     = _z.size();
      func_ineq_data.s        = _s;

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
          result[i] = -coeff[i] + std::inner_product(x, x+n-1, bmat.begin() + i*n, 0.0);
        }
      };

      auto n_ineq = _z.size();
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
        assert(m == data.npts+s+1);
        assert(n == s+2);
        using std::conj;
        using std::real;
        if (grad)
        {
          for (size_t i = 0; i < s+1; i++)
          {
            for (size_t j = 0; j < s+2; j++)
              grad[i*n+j] = 0;
            grad[i*n+i] = -1;
            result[i] = -x[i];
          }
          for (size_t i = 0; i < data.npts; i++)
          {
            const auto g = 
              std::inner_product(x, x+n-1, cmat.begin() + i*n,complex_type{0});
            for (size_t j = 0; j < s+1; j++)
            {
              const auto df = cmat[i*n+j]*conj(g);
              grad[(s+1+i)*n+j] = real(df + conj(df));
            }
            grad[(s+1+i)*n+s+1] = -1;

            const auto re = real(g*conj(g));
            result[s+1+i] = (re-1)-x[n-1];
          }
        }
        else
        {
          for (size_t i = 0; i < s+1; i++)
            result[i] = -x[i];

          for (size_t i = 0; i < data.npts; i++)
          {
            const auto g = 
              std::inner_product(x, x+n-1, cmat.begin() + i*n,complex_type{0});
            const auto re = real(g*conj(g));
            result[s+1+i] = (re-1)-x[n-1];
          }
        }
      };
  
      auto method = nlopt::LD_SLSQP;
      nlopt::opt opt(method, _s+2);

      opt.set_min_objective(func,NULL);

      const real_type tol = _tol;
      opt.add_equality_mconstraint(
          func_eq, 
          &func_eq_data, 
          std::vector<real_type>(_p+1,tol));

      std::vector<real_type> &x = _solution;
      x.resize(_s+2, 1);

      if (x_guess.empty())
      {
        std::generate(x.begin(), x.end(), [](){ return drand48(); });
//        std::fill(x.begin(), x.end(), 1);
      }
      else
      {
        std::copy(x_guess.begin(), x_guess.end(), x.begin());
        x.back() = 0;
      }
      opt.add_inequality_mconstraint(
          func_ineq, 
          &func_ineq_data, 
          std::vector<real_type>(_s+1+n_ineq,tol)); 

//      opt.set_xtol_rel(tol);
      opt.set_xtol_abs(tol);
//     opt.set_ftol_abs(tol);
//      opt.set_ftol_rel(tol);

      real_type& minf = _fmin;
      decltype(opt.optimize(x,minf)) result;
      try
      {
        result = opt.optimize(x,minf);
      }
      catch (std::exception &e)
      {
        if (_verbose)
          printf(std::cerr, " Exception caught: % \n", e.what());
      }

      if (_verbose)
      {
//        printf(std::cerr, "result= % \n", result);
        printf(std::cerr, "minf= % \n", minf);
      }

      if (_verbose)
      {
        const auto &poly = x;
        const auto &bmat = _bmat;
#if 0
        std::cerr << "fixed_coeff_comp= " << std::setprecision(16);
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          std::cerr << res << " ";
        }
        std::cerr << std::endl;

        std::cerr << "fixed_coeff_gold= " ;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = fixed_coeff[i];
          std::cerr << res << " ";
        }
        std::cerr << std::endl;
#endif
        std::cerr << "fixed_coeff_diff= " ;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          std::cerr << std::abs(res-fixed_coeff[i])/fixed_coeff[i] << " ";
        }
        std::cerr << std::endl;
//        std::cerr << " ===================================== \n";
      }
      real_type norm = 0;
      {
        const auto &poly = x;
        const auto &bmat = _bmat;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          const auto f = (res-fixed_coeff[i])/fixed_coeff[i];
          norm += f*f;
        }
        norm = std::sqrt(norm/(_p+1));
      }
      return norm;

    }
};

template<typename real_type>
static auto linspace(const real_type a, const real_type b, const size_t n)
{
  std::vector<real_type> res(n);
  std::iota(res.begin(), res.end(), 0);
  for (auto& x : res)
    x = a + (b-a)*x/(n-1); 
  return res;
}


template<typename real_type, typename complex_type>
static auto maximizeH(const size_t p, const size_t s,const std::vector<complex_type>& ev_space)
{
#if 0
  using Optimizer = OptimizerT<real_type,complex_type>;
#else
  using Optimizer = OptimizerL_CMP<real_type,complex_type>;
#endif

  Optimizer opt(s,p,ev_space);

  real_type h_min = 0; 
  real_type h_max = 2.01*s*s*std::abs(opt.get_zmin())/p;

  const auto max_steps = 1000;
  const auto tol_bisect = 0.01;

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

    opt.unset_verbose();
    opt.optimize(h);
    printf(std::cerr, "step= %  h_min= %  h_max= %  -- h= %  val= % \n",
        step, h_min, h_max, h, opt.get_fmin());

    if (std::abs(opt.get_fmin()) < 1.0e-12) //std::abs(get<1>(res)) < 1.0e-12)
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
  opt.set_verbose();
  opt.optimize(h);
  printf(std::cerr, "solution:  h_min= %  h_max= %  -- h= %  val= % \n",
       h_min, h_max, h, opt.get_fmin());

  return std::make_tuple(opt.get_solution(), h, opt);
}

static auto maximizeHdriver(int order, int stages, int npoints)
{
  using std::get;
  using namespace std::literals;
  using real_type    =  double;
  using complex_type = std::complex<real_type>;

  auto p = order;
  auto npts = npoints;

  auto beta   = 1.00;
  auto s = stages;

  /*  build eigen-value space */

  const auto imag_lim = std::abs(beta);
  const auto l1 = linspace(0.0,imag_lim,npts/2);
  const auto l2 = linspace(-1.0,0.0,npts);

  std::vector<complex_type> ev_space;
  if (imag_lim > 0)
    for (auto& x : l1)
      ev_space.emplace_back(1i*x);
  for (auto& x : l2)
    ev_space.emplace_back(1i*imag_lim + x);
  if (imag_lim > 0)
    for (auto &x : l1)
      ev_space.emplace_back(-1.0 + 1i*x);
  
  printf(std::cerr, "p= % \n", p);
  printf(std::cerr, "s= % \n", s);
  printf(std::cerr, "npts= % \n", npts);


  auto res = maximizeH<real_type,complex_type>(p,s,ev_space);


  const auto & h    = get<1>(res);
#if 0
  const auto & poly = get<0>(res);
  std::cout << "coeff = { \n";
  for (size_t i = 0; i < poly.size() -1 ; i++)
  {
    std::cout << std::setprecision(16) << poly[i] << ", ";
  }
  std::cout << poly.back()  << " };\n";
#endif
  std::cerr << "h= " << h << std::endl;
  std::cerr << "h/s^2= " << h/(s*s) << std::endl;
//  std::cout << "h_imag= " << beta << std::endl;

  return std::make_tuple(h,std::move(get<2>(res)));
}

template<typename Optimizer>
static auto minimizeS(
    Optimizer &opt, 
    const size_t smax,
    const double h_real, const double h_imag)
{
  auto s_min = opt.get_p() + 1;
  auto s_max = smax;

  auto converged = false;
  auto s = s_max;
  while (1)
  {
    if (s_max - s_min <= 1)
    {
      if (converged)
        break;
      else
        s = s_max;
    }
    else
    {
      if (s_max - s_min > 2)
        s = std::max(s_min,static_cast<decltype(s)>(0.25*s_min + 0.75*s_max+0.5));
      else
        s = (s_min + s_max)>>1;
    }
    opt.set_s(s);
    
    opt.unset_verbose();
    opt.optimize(h_real, h_imag);
//    printf(std::cerr, " s_min= %  s_max= %  -- s= %  val= % \n", s_min, s_max, s, opt.get_fmin());

    if (std::abs(opt.get_fmin()) < 1.0e-12) 
    {
      converged = true;
      s_max = s;
    }
    else
    {
      converged = false;
      s_min = s;
    }
  }

  s = std::min(smax, static_cast<decltype(s)>(s*1.05));
  assert(converged);
  return s;

  opt.set_verbose();
  opt.set_s(s);
  opt.optimize(h_real, h_imag);
  printf(std::cerr, " solution: s_min= %  s_max= %  -- s= %  val= % \n",  s_min, s_max, s, opt.get_fmin());

  return s;
}

template<size_t ORDER>
static auto orderN(
    const int stages, const int npts,
    const double nodes_[ORDER],
    const double nodesC_[ORDER-1])
{
  std::vector<std::vector<double>> coeff;
  using std::get; 


  auto res= maximizeHdriver(ORDER,stages,npts);
  const auto h_base = get<0>(res);
  auto& opt = get<1>(res);
  const auto s_base = opt.get_s();

  std::vector<double> nodes, nodesC;
  for (int i = 0; i < ORDER; i++)
  {
    auto x= nodes_[i];
    x = (1+x)/2;
    nodes.push_back(x);
  }
  for (int i = 0; i < ORDER-1; i++)
  {
    auto x= nodesC_[i];
    x = (1+x)/2;
    nodesC.push_back(x);
  }
 
  const auto smax = opt.get_s();
  for (int i = 0; i < ORDER; i++)
  {
    auto h = h_base*nodes[i];
    auto s = minimizeS(opt, smax, h, nodes[i]);
//    auto s = smax;
    opt.set_s(s);
    opt.unset_verbose();
    const auto norm = opt.optimize(h, nodes[i]);
    printf(std::cerr, " %: i= %  s= % node= %  fmin= %  norm= %\n", ORDER, i, s, nodes[i], opt.get_fmin(), norm);
    const auto & poly = opt.get_solution();
    coeff.push_back(poly);
  }
  std::cerr << " -- Embedded method -- " << std::endl;
  opt.set_p(ORDER-1);
  for (int i = 0; i < ORDER-1; i++)
  {
    auto h = h_base*nodesC[i];
    auto s = minimizeS(opt, smax, h, nodesC[i]);
//    auto s = smax;
    opt.set_s(s);
    opt.unset_verbose();
    const auto norm = opt.optimize(h, nodesC[i]);
    printf(std::cerr, " %: i= %  s= % node= %  fmin= %  norm= %\n", ORDER-1,i, s, nodesC[i], opt.get_fmin(), norm);
    const auto & poly = opt.get_solution();
    coeff.push_back(poly);
  }
  printf(std::cerr, "---------------------------\n");

  return std::make_tuple(h_base, coeff);
}


static auto order8(const int stages, const int npts)
{
  assert(stages >= 8);
  const double nodes [8] = {-0.960289856497536231684, -0.7966664774136267395916, -0.5255324099163289858177, -0.1834346424956498049395, 0.1834346424956498049395, 0.525532409916328985818, 0.796666477413626739592, 0.9602898564975362316836};
  const double nodesC[7] = {-0.9491079123427585245262,-0.7415311855993944398639,-0.4058451513773971669066,0,0.4058451513773971669066,0.7415311855993944398639,0.9491079123427585245262};

  return orderN<8>(stages, npts, nodes,nodesC);
}

static auto order4(const int stages, const int npts)
{
  assert(stages >= 4);
  const double nodes [4] = {-0.8611363115940525752239, -0.3399810435848562648027, 0.3399810435848562648027, 0.8611363115940525752239}; 
  const double nodesC[3] = {-0.7745966692414833770359, 0, 0.7745966692414833770359};
  return orderN<4>(stages, npts, nodes,nodesC);
}

static auto order2(const int stages, const int npts)
{
  assert(stages >= 2);
  const double nodes [2] = {-0.5773502691896257645091, 0.5773502691896257645091};
  const double nodesC[1] = {0};
  return orderN<2>(stages, npts, nodes,nodesC);
}

template<typename Func>
static auto compute_coeff(Func func, int stage_min, int stage_max)
{
  const auto npts = 1000;
  using coeff_type = decltype(func(stage_max, npts));

  std::vector<coeff_type> coeff;

  for (int s = stage_min; s <= stage_max; s++)
    coeff.push_back(func(s,npts));

  return coeff;

}

void run_loop(const int order, const int min_stages, const int max_stages)
{
  printf(std::cerr, " -------------- \n");
  printf(std::cerr, " order= %  min_stages= %  max_stages= %\n", order, min_stages, max_stages);
  printf(std::cerr, " -------------- \n");

  assert(order == 4 || order == 8 || order == 2);
  auto func = order4;
  if (order == 2)
    func = order2;
  else if  (order == 8)
    func = order8;

  auto res = compute_coeff(func, min_stages, max_stages);
 
 
  using std::cout;
  using std::endl;
  using std::get;
  cout << std::setprecision(16);
  cout << "const int order = " << order << ";\n";
  cout << "const int nelements= " << res.size() << ";\n";
  cout << "const double h_base = {\n";
  for (const auto& x : res)
  {
    const auto h = get<0>(x);
    cout << h << ", " << endl;
  }
  cout << "0}; " << endl;
#if 1
  cout << "const int offset = { \n";
  auto sum = 0;
  for (const auto& x : res)
  {
    cout << sum << ", " << endl;
    const auto c = get<1>(x);
    for (const auto &cc : c)
      sum += cc.size();
  }
  cout << sum << "}; " << endl;
#endif
  cout << "const int sizes = { \n";
  for (const auto& x : res)
  {
    const auto c = get<1>(x);
    cout << " ";
#if 1
    for (size_t i = 0; i < c.size()-1; i++)
      cout << c[i].size() << ", ";
#endif
    cout << c.back().size() << ",\n";
  }
  cout << "0}; " << endl;
  cout << "const double coeff = { \n";
  for (const auto& x : res)
  {
    const auto c = get<1>(x);
    for (auto &cc : c)
      for (auto &y : cc)
      {
        if (std::abs(y) < 1.0e-14)
          cout << 0 << ", ";
        else
          cout << y << ", ";
      }
    cout << endl;
  }
  cout << "0}; " << endl;
}

void single_order(const int order, int min_stages, int max_stages)
{
  printf(std::cerr, " -------------- \n");
  printf(std::cerr, " order= %  min_stages= %  max_stages= %\n", order, min_stages, max_stages);
  printf(std::cerr, " -------------- \n");
  
  using std::get;
  using vector = std::vector<double>;
  std::vector<double> h_base;
  std::vector<std::tuple<vector,vector>> coeff;

  int npts = 1000;
  for (int s = min_stages; s <= max_stages; s++)
  {
    printf(std::cerr, " >>>>>> s= % <<<<<<<<< \n", s);
    auto res = maximizeHdriver(order,s,npts);
    const auto h = get<0>(res);
    auto& opt = get<1>(res);
    const auto norm1 = opt.optimize(h);
    const auto fmin1 = opt.get_fmin();
    h_base.emplace_back(h);
    auto poly1 = opt.get_solution();
    opt.set_p(order-1);
    const auto norm2 = opt.optimize(h);
    const auto fmin2 = opt.get_fmin();
    printf(std::cerr,  "s= %: fmin1= %  fmin2= %   norm1= %  norm2= %\n", order, fmin1, fmin2, norm1, norm2);
    auto poly2 = opt.get_solution();
    coeff.emplace_back(poly1, poly2);
  }
  assert(coeff.size() == h_base.size());
  using std::cout;
  using std::endl;
  using std::get;
  cout << std::setprecision(16);
  cout << "const int order = " << order << ";\n";
  cout << "const int nelements= " << coeff.size() << ";\n";
  cout << "const double h_base = {\n";
  for (const auto& h : h_base)
    cout << h << ", " << endl;
  cout << "0}; " << endl;
  cout << "const int offset = { \n";

  auto sum = 0;
  for (const auto& c : coeff)
  {
    cout << sum << ", " << endl;
    assert(get<0>(c).size() == get<1>(c).size());
    sum += get<0>(c).size()*2;
  }
  cout << sum << "}; " << endl;
  cout << "const double coeff = { \n";
  for (const auto& c : coeff)
  {
    auto outx = [](const auto x) 
    {
      if (std::abs(x) < 1.0e-14)
        cout << 0 << ", ";
      else
        cout << x << ", ";
    };
    for (auto &x: get<0>(c))
      outx(x);
    for (auto &x: get<1>(c))
      outx(x);
    cout << endl;
  }
  cout << "0}; " << endl;
}

int main(int argc, char * argv[])
{
  srand48((size_t)argv+1234);
  using std::get;
#if 1
  assert(argc >= 2);
  const auto min_stages = atoi(argv[1]);
  const auto max_stages = atoi(argv[2]);
  single_order(4, min_stages, max_stages);
#elif 1
  assert(argc >= 2);
  const auto min_stages = atoi(argv[1]);
  const auto max_stages = atoi(argv[2]);
  run_loop(4, min_stages, max_stages);
#else
  auto order = 4;
  auto stages = 100;
  auto npts = 1000;
  auto res = maximizeHdriver(order,stages,npts);
  auto h = get<0>(res);
  auto &opt = get<1>(res);
  const auto &poly = opt.get_solution();
  std::cout << "coeff = { \n";
  auto sum = 0.0;
  for (size_t i = 0; i < poly.size() -1 ; i++)
  {
    auto x = poly[i];
    if (std::abs(x) < 1.0e-14) x = 0;
    std::cout << std::setprecision(16) << x << ", ";
    sum += x;
  }
  {
    auto x = poly.back();
    if (std::abs(x) < 1.0e-14) x = 0;
    sum += x;
    std::cout << x  << " };\n";
  }
  std::cerr << "h= " << h << std::endl;
  std::cerr << "h/s^2= " << h/(stages*stages) << std::endl;
  std::cerr << "sum= " << sum << std::endl;
#endif
  return 0;
}
