#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <cassert>
#include "common.h"

template<size_t N, typename T, typename Real>
class ExpansionBaseT
{
  protected:
    using vector_type = std::array<Real,N>;
    using matrix_type = std::array<vector_type,N>;
  public:
    using storage     = std::array<T,N>;
    static constexpr auto size() { return N; }
};

template<size_t N, typename T, typename Real>
class ExpansionT;

template<typename T, typename Real>
class ExpansionT<1,T,Real> : public ExpansionBaseT<1,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<1,T,Real>;
    static constexpr size_t N = 1;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    {
      return Real{1};
    }

    static constexpr auto weight(const size_t i) 
    {
      return Real{1};
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    { 
      return Real{0.5};
    }
    static constexpr auto maxAbsMEV() 
    {
      return Real{2};
    }
    static constexpr auto maxAbsPEV() 
    { 
      return Real{2};
    }
};

template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  protected:
    using base_type = ExpansionBaseT<2,T,Real>;
    static constexpr size_t N = 2;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {3.0,  0.4641016151377546},
        {-6.464101615137754,  3.0}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] = {0.5,0.5};
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.250000000000000000000,  0},
        {0.53867513459481288225, 0.250000000000000000000}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV = 3.5;
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV = 0.25;
      return maxAbsPEV; 
    }
};

template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  protected:
    using base_type  = ExpansionBaseT<3,T,Real>;
    static constexpr size_t N = 3;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {5.0,  1.1639777949432226, -0.1639777949432225},
        {-5.727486121839514,  2.0, 0.7274861218395141},
        {10.163977794943223,  -9.163977794943223, 5.0}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.2777777777777778,
        0.4444444444444444,
        0.2777777777777778
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.13888888888888888889, 0, 0},
        {0.30026319498086459244, 0.22222222222222222222,  0},
        {0.26798833376246945173, 0.48042111196938334790,  0.13888888888888888889}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{5.1};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.22222};
      return maxAbsPEV; 
    }
};


template<size_t ORDER, typename PDE>
struct ODESolverT
{
  using Real           = typename PDE::Real;
  using Vector         = typename PDE::Vector;
  using Expansion      = ExpansionT<ORDER,Vector,Real>;
  using range_iterator = make_range_iterator<size_t>;


  PDE _pde;
  Real _time;
  typename Expansion::storage _x, _rhs;

  auto expansionRange() const 
  {
    return make_range_iteratorT<0,Expansion::size()>{};
  }

  auto time() const {return _time;}

  ODESolverT(const PDE &pde) : _pde{pde}, _time{0}
  {
    for (auto k : expansionRange())
    {
      _x  [k].resize(_pde.resolution());
      _rhs[k].resize(_pde.resolution());
    }
  };

  PDE& pde() { return _pde; }
  const PDE& pde() const { return _pde; }


  void rhs(const Vector &u0)
  {
    using std::get;
    const auto& x = _x;
    auto& rhs = _rhs;

    /* compute RHS */
    for (auto k : expansionRange())
    {
      _pde.compute_rhs(rhs[k], x[k]);

      assert(_x[k].size() == u0.size());
      for (auto l : expansionRange())
        for (auto v : make_zip_iterator(rhs[k], x[l], u0))
        {
          get<0>(v) += Expansion::matrix(k,l) * (get<2>(v) - get<1>(v));
        }
    }

    /* precondition RHS */
    for (auto i : range_iterator{0,u0.size()})
    {
      std::array<Real,Expansion::size()> tmp;
      for (auto& t : tmp)
        t = 0;
        
      for (auto k : expansionRange())
        for (auto l : expansionRange())
        {
          tmp[k] += Expansion::preconditioner(k,l)*_rhs[l][i];
        }

      for (auto k : expansionRange())
        _rhs[k][i] = tmp[k];
    }
  }
  
  void iterate(const Vector &u0)
  {
    rhs(u0);

    const Real omega = Real{0.7}/(Expansion::maxAbsPEV()*(Expansion::maxAbsMEV() + _pde.AbsEV()));
#if 0
    printf(std::cerr, "omega= %   PEV= %  MEV= %  cfl= % \n",
        omega,
        Expansion::maxAbsPEV(),
        Expansion::maxAbsMEV(),
        _pde.cfl());
#endif

    /* use LegendreP(1,z) for update */
    using std::get;
    for (auto k : expansionRange())
      for (auto v : make_zip_iterator(_x[k], _rhs[k]))
      {
        get<0>(v) = get<0>(v) + 2.0*omega*get<1>(v);
      }
  }

  void solve_system(const Vector& u0)
  {
    using std::get;
    constexpr auto niter = 40;
    std::array<Real,Expansion::size()> error;
    for (auto iter : range_iterator{0,niter})
    {
      for (auto k : expansionRange())
      {
        _pde.apply_bc(_x[k]);
        error[k] = 0;
      }

      iterate(u0);
      
      Real err = 0;
      int cnt = 0;
      constexpr auto eps = Real{1.0e-7};
      for (auto k : expansionRange())
      {
        for (auto i : range_iterator{1,u0.size()-1})
        {
          error[k] += square(_rhs[k][i]); ///(_x[k][i] + eps));
          cnt += 1;
        }
        err += std::max(err,std::sqrt(error[k]/cnt));
      }
      if (iter == niter - 1)
        printf(std::cerr, "   >> iter= %  err= % \n", iter, err);
    }
  }

  void update()
  {
    using std::get;
    static auto du = _pde.state();

    for (auto k : expansionRange())
    {
      _x[k] = _pde.state();
    }
    solve_system(_pde.state());

    for (auto k : expansionRange())
    {
      _pde.apply_bc(_x[k]);
      _pde.compute_rhs(_rhs[k], _x[k]);
    }

    std::fill(du.begin(), du.end(), 0);
    for (auto k : expansionRange())
    {
      const auto w = Expansion::weight(k);
      for (auto v : make_zip_iterator(du,_rhs[k]))
      {
        get<0>(v) += w*get<1>(v);
      }
    }

#if 0
    for (auto u : du)
      if (std::abs(u) > 1.0e-10)
        printf(std::cerr,  "u= % ", u);
#endif
    _pde.update(du);
    _time += _pde.dt();
  }

#if 0
  void integrate(Real dt)
  {
  }
#endif

};

template<typename real_type>
struct PDEDiffusion
{
  using Real   = real_type;
  using Vector = std::vector<Real>;
  using range_iterator = make_range_iterator<size_t>;

  Vector _f;

  Real _cfl;
  Real _dx;
  Real _diff;
  Real _zeta;

  Real dt() const {return _cfl * 0.5*square(_dx)/_diff;}
  Real AbsEV() const
  {
    return dt() * 4.0*_diff/square(_dx);  /* 2.0 * cfl */
  }

  auto dx() const { return _dx; }

  PDEDiffusion(const size_t n) : _f{Vector(n)}
  {
  }

  static void periodic_bc(Vector &f) 
  {
    const auto n = f.size();
    f[0  ] = f[n-2];
    f[n-1] = f[1  ];
  }

  static void free_bc(Vector &f) 
  {
    const auto n = f.size();
    f[0  ] = f[  1];
    f[n-1] = f[n-2];
  }

  auto cfl() const { return _cfl; }
  auto resolution() const { return _f.size(); }

  void apply_bc(Vector &f) const
  {
//    periodic_bc(f);
    free_bc(f);
  }

  const Vector& state() const { return _f; }

  void update(const Vector &df) 
  {
    using std::get;
    for (auto v: make_zip_iterator(_f,df))
    {
      get<0>(v) += get<1>(v);
    }
  }

  template<typename Func>
  void compute_rhs(Vector &res, const Vector &x, Func func)
  {
    const auto c = dt() * _diff/square(_dx);
    for (auto i : range_iterator{1,x.size() - 1})
    {
      res[i] = c * (x[i+1] - Real{2.0} * x[i] + x[i-1]);
      res[i] = func(res[i]);
    }
  }
  void compute_rhs(Vector &res, const Vector &x)
  {
    compute_rhs(res, x, [](const auto x) { return x; });
  }

  void set_ic()
  {
    using std::max;

    /* set  a profile of delta function */

    auto &f = _f;

    const int n = f.size();
    const auto dx = _dx;
    const auto L = dx*(n-2);
    const auto ic = n>>1;

    const auto ampl = Real{1.0};


    const int di = 10;
    std::fill(f.begin(), f.end(), 0);
    for (int i = -di; i <= di; i++)
    {
      const auto fi = max(ampl*(di - abs(i)),Real{0});
      f[ic + i] = fi;
    }
    std::fill(f.begin(), f.end(), 0);
    for (int i = -di; i <= di; i++)
    {
      const auto fi = max(ampl*(di - abs(i)),Real{0});
      f[ic - ic/2 + i] = fi;
    }
    for (int i = -di; i <= di; i++)
    {
      const auto fi = max(ampl*(di - abs(i)),Real{0});
      f[ic + ic/3 + i] = fi;
    }

  }
};

template<typename Solver>
void dump2file(const Solver &solver, std::string fileName = std::string{})
{
  const auto &f = solver.pde().state();
  const auto n = f.size();
  const auto dx = solver.pde().dx();
  std::ofstream foutFn(fileName);
  std::ostream &fout = fileName.empty() ? std::cout : foutFn;
  fout << "# time= " << solver.time() << std::endl;
  fout << "# n= " << n-2 << std::endl;
  for (size_t i = 1; i < n-1; i++)
  {
    const auto x = (i-1+0.5)*dx/(n-2);
    fout << x << " " << std::setprecision(16) << f[i] << std::endl;
  }
}

int main(int argc, char * argv[])
{
  using Real = double;

  const size_t ncell = argc > 1 ? atoi(argv[1])+2 : 128+2;
  printf(std::cerr, "ncell= %\n", ncell);

  const Real tau = argc > 2 ? atof(argv[2]) : 1.0;
  printf(std::cerr, "tau= %\n", tau);
  


  constexpr auto ORDER = 3;
  using PDE = PDEDiffusion<Real>;
  using Solver = ODESolverT<ORDER,PDE>;


  Solver solver(PDE{ncell});

  solver.pde()._dx   = 1;
  solver.pde()._diff = 1;
  solver.pde()._cfl  = 0.5*8;  /* stable for cfl <= 0.5 */

  const auto dt = solver.pde().dt();
  const size_t nstep = std::max(size_t{1}, static_cast<size_t>(tau/dt));

  printf(std::cerr, "dt= %  tau= %   nstep= %\n",
      dt, tau, static_cast<int>(tau/dt));

  solver.pde().set_ic();
  dump2file(solver, "ic.txt");
  for (size_t step = 1; step <= nstep; step++)
  {
    solver.update();
    printf(std::cerr, "step= % : time= % \n", step, solver.time());
  }
  dump2file(solver);




  return 0;  

}




