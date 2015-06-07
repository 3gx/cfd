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

    void optimize(const real_type h_real , const real_type h_imag = 1, const std::vector<real_type> x_guess = std::vector<real_type>())
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
        std::fill(x.begin(), x.end(), 1);
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

      opt.set_xtol_abs(tol);

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
        printf(std::cerr, "result= % \n", result);
        printf(std::cerr, "minf= % \n", minf);
      }

      if (_verbose)
      {
        const auto &poly = x;
        const auto &bmat = _bmat;
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

        std::cerr << "fixed_coeff_diff= " ;
        for (int i = 0; i < _p+1; i++)
        {
          const auto res = std::inner_product(poly.begin(), poly.end()-1, bmat.begin() + i*(_s+2), 0.0);
          std::cerr << std::abs(res-fixed_coeff[i])/fixed_coeff[i] << " ";
        }
        std::cerr << std::endl;
        std::cerr << " ===================================== \n";
      }

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
  using Optimizer = OptimizerT<real_type,complex_type>;

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

    if (opt.get_fmin() < 1.0e-12) //std::abs(get<1>(res)) < 1.0e-12)
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

static auto maximizeHdriver(int order, double stiffness, int npoints, double beta_scale = 1.0)
{
  using std::get;
  using namespace std::literals;
  using real_type    =  double;
  using complex_type = std::complex<real_type>;

  auto p = order;
  auto npts = npoints;
  auto stiff_kappa = stiffness;

  auto beta   = beta_scale/0.95; 
  auto alpha_s = 0.15;
  auto s = std::max(p,static_cast<int>(std::sqrt(stiff_kappa/alpha_s) + 1));
#if 0
  if (s < 12)
  {
    alpha_s = 0.07;
    s = std::max(p,static_cast<int>(std::sqrt(stiff_kappa/alpha_s) + 1));
  }
#endif

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
  printf(std::cerr, "kappa_stiff= % \n", stiff_kappa);
  printf(std::cerr, "beta= % \n",  beta);


  auto res = maximizeH<real_type,complex_type>(p,s,ev_space);


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
  std::cout << "h_imag= " << beta << std::endl;

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
    printf(std::cerr, " s_min= %  s_max= %  -- s= %  val= % \n",
        s_min, s_max, s, opt.get_fmin());

    if (opt.get_fmin() < 1.0e-12) 
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
  opt.set_verbose();
  opt.set_s(s);
  opt.optimize(h_real, h_imag);
  printf(std::cerr, " solution: s_min= %  s_max= %  -- s= %  val= % \n",
       s_min, s_max, s, opt.get_fmin());

  return s;
}


static void order8(const int stiffness, const int npts)
{
  constexpr int order = 8;
  auto node = [](const int i) 
  {
    static const double nodes[order] = {
      -0.960289856497536231684,
      -0.7966664774136267395916,
      -0.5255324099163289858177,
      -0.1834346424956498049395,
      0.1834346424956498049395,
      0.525532409916328985818,
      0.796666477413626739592,
      0.9602898564975362316836
    };
    return  (nodes[i]+1)*0.5;
  };

  using std::get;
  auto res = maximizeHdriver(order,stiffness,npts);
  const auto& h_max  = get<0>(res);
  auto& opt = get<1>(res);

  const auto h_step  = h_max*0.95;
  opt.set_verbose();
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "h_step= " << h_step << std::endl;
  for (int i = 0; i < order; i++)
  {
    auto h_try = h_step * node(i);
    assert(h_try > 0);
    assert(h_try < h_step);
    opt.unset_verbose();
    auto node_im = node(i)*1.0;
    opt.optimize(h_try, node_im);
    opt.set_verbose();
    opt.optimize(h_try,node_im,opt.get_solution());
    std::cout << "h= " << h_try << std::endl;
    std::cout << "coeff[" << i<< "] =\n{ \n";
    const auto & poly = opt.get_solution();
    const auto & h    = h_try;
    for (size_t i = 0; i < poly.size() -1 ; i++)
    {
      std::cout << std::setprecision(16) << poly[i] << ", ";
    }
    std::cout << poly.back()  << " };\n\n";
//    std::cout << "h/s^2= " << h_try/(opt.get_s()*opt.get_s()) << std::endl;
  }
}

#if 0
static void order8_1(const int stiffness, const int npts)
{
  constexpr int order = 8;
  const double nodes_[order] = {
    -0.960289856497536231684,
    -0.7966664774136267395916,
    -0.5255324099163289858177,
    -0.1834346424956498049395,
    0.1834346424956498049395,
    0.525532409916328985818,
    0.796666477413626739592,
    0.9602898564975362316836
  };
  std::vector<double> nodes;
  for (auto x : nodes_)
    nodes.push_back((1+x)/2);
  
  const double stiffy[order] = 
  {
    2*nodes[0],
    2*nodes[1],
    2*nodes[2],
    2*nodes[3],
    2*nodes[4],
    nodes[5],
    nodes[6],
    nodes[7]
  };

   using std::get; 
  const auto res_base  = maximizeHdriver(order,stiffness,npts);
  const auto  h_base   = get<0>(res_base);

  for (int i = 0; i < order; i++)
  {
    using std::get;
    const auto kappa = stiffness* std::min(1.0,stiffy[i]);
    auto res = maximizeHdriver(order,kappa,npts,nodes[i]);
    auto   h1  = get<0>(res);
    auto   h  = h_base*nodes[i];
#if 0
    int  dstiff =1;
    while (h > h1)
    {
      dstiff++;
      res = maximizeHdriver(order,stiffness*nodes[i]+dstiff,npts,nodes[i]);
      h1 = get<0>(res);
      auto& opt = get<1>(res);
      printf(std::cout, "step= %  node= % h= % h1= % s= % \n", i, nodes[i], h, h1, opt.get_s());
    }
#endif
    auto& opt = get<1>(res);
    printf(std::cout, "------------------------\n");
    printf(std::cout, "step= %  node= % h= % h1= % s= % \n", i, nodes[i], h, h1, opt.get_s());
    assert(h < h1);
    opt.optimize(h,1,opt.get_solution());

    std::cout << "h= " << h<< std::endl;
    std::cout << "coeff[" << i<< "] =\n{ \n";
    const auto & poly = opt.get_solution();
    for (size_t i = 0; i < poly.size() -1 ; i++)
    {
      std::cout << std::setprecision(16) << poly[i] << ", ";
    }
    std::cout << poly.back()  << " };\n\n";
  }
}
#endif

static auto order8_2(const int stiffness, const int npts)
{
  std::vector<std::vector<double>> coeff;
  using std::get; 

  constexpr int order = 8;


  auto res= maximizeHdriver(order,stiffness,npts);
  const auto h_base = get<0>(res);
  auto& opt = get<1>(res);
  const auto s_base = opt.get_s();

  const double nodes8_[8] = {-0.960289856497536231684, -0.7966664774136267395916, -0.5255324099163289858177, -0.1834346424956498049395, 0.1834346424956498049395, 0.525532409916328985818, 0.796666477413626739592, 0.9602898564975362316836};
  const double nodes7_[7] = {-0.9491079123427585245262,-0.7415311855993944398639,-0.4058451513773971669066,0,0.4058451513773971669066,0.7415311855993944398639,0.9491079123427585245262};
  std::vector<double> nodes8, nodes7;
  for (auto x : nodes8_)
  {
    x = (1+x)/2;
    nodes8.push_back(x);
  }
  for (auto x : nodes7_)
  {
    x = (1+x)/2;
    nodes7.push_back(x);
  }
 
  const auto smax = opt.get_s();
  for (int i = 0; i < order; i++)
  {
    auto s = minimizeS(opt, smax, h_base*nodes8[i], nodes8[i]);
    auto h = h_base*nodes8[i];
    opt.set_s(s);
    opt.optimize(h, nodes8[i]);
    printf(std::cout, "------------------------\n");
    std::cerr << "s= " << s << " node= " << nodes8[i] <<  std::endl;
    std::cout << " h= " << h<<  std::endl;
    std::cout << "coeff[" << i<< "] =\n{ \n";
    const auto & poly = opt.get_solution();
    for (size_t i = 0; i < poly.size()-1; i++)
    {
      std::cout << std::setprecision(16) << poly[i] << ", ";
    }
    std::cout << poly.back()  << " };\n\n";
  }
}

void order4()
{
  const double nodes4_[4] = {-0.8611363115940525752239, -0.3399810435848562648027, 0.3399810435848562648027,0.8611363115940525752239};
  const double nodes3_[3] = {-0.7745966692414833770359, 0, 0.7745966692414833770359};
  std::vector<double> nodes4, nodes3;
  for (auto x : nodes4_)
  {
    x = (1+x)/2;
    nodes4.push_back(x);
  }
  for (auto x : nodes3_)
  {
    x = (1+x)/2;
    nodes3.push_back(x);
  }
 
}

void order2()
{
  auto node = [](const int i) 
  {
    const double nodes[2] = {-0.5773502691896257645091, 0.5773502691896257645091};
    return  (nodes[i]+1)*0.5;
  };
}

int main(int argc, char * argv[])
{
  const auto npts = 1000;
  const auto stiffness = 100;
#if 1
  order8_2(stiffness,npts);
#else
  const auto order = 8;

  using std::get;
  auto res = maximizeHdriver(order,stiffness,npts);

  auto& h_max = get<0>(res);
  auto& opt   = get<1>(res);

  printf(std::cerr, " ------------ \n" );

  {
    auto h_try = floor(h_max*0.95);
    opt.set_verbose();
    opt.optimize(h_try);

    std::cout << "h= " << h_try << std::endl;
    std::cout << "coeff = { \n";
    const auto & poly = opt.get_solution();
    const auto & h    = h_try;
    for (size_t i = 0; i < poly.size() -1 ; i++)
    {
      std::cout << std::setprecision(16) << poly[i] << ", ";
    }
    std::cout << poly.back()  << " };\n\n";
  }
#endif
  return 0;
}
