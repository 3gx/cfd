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
    using vector_type = typename base_type::vector_type;
    using matrix_type = typename base_type::matrix_type;


  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    {
      const matrix_type _matrix{vector_type{2}};
      return _matrix[j][i]; 
    }

    static constexpr auto weight(const size_t i) 
    {
      const vector_type _weight{1};
      return _weight[i];  
    }
    static constexpr auto zero(const size_t i) 
    { 
      const vector_type _zero{1};
      return _zero[i]; 
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    { 
      const matrix_type _preconditioner{vector_type{0.5}};
      return _preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real _maxAbsMEV{2};
      return _maxAbsMEV; 
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real _maxAbsPEV{0.5};
      return _maxAbsPEV;
    }
};

template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<2,T,Real>;
    using vector_type = typename base_type::vector_type;
    using matrix_type = typename base_type::matrix_type;

    static constexpr matrix_type _matrix{
      vector_type{3.0,  0.4641016151377546},
      vector_type{-6.464101615137754,  3.0}
    };

    static constexpr vector_type _weight{Real{0.5},Real{0.5}};
    static constexpr vector_type _zero{1,1};
    
    static constexpr matrix_type _preconditioner{
      vector_type{0.250000000000000000000,  0},
      vector_type{0.53867513459481288225, 0.250000000000000000000}
    };

    static constexpr Real _maxAbsMEV{3.5};
    static constexpr Real _maxAbsPEV{0.25};

  public:
    static constexpr auto matrix(const size_t i, const size_t j) { return _matrix[j][i]; }
    static constexpr auto weight(const size_t i)  { return _weight[i];  }
    static constexpr auto zero(const size_t i) { return _zero[i]; }
    static constexpr auto preconditioner(const size_t i, const size_t j) { return _preconditioner[j][i]; }
    static constexpr auto maxAbsMEV() { return _maxAbsMEV; }
    static constexpr auto maxAbsPEV() { return _maxAbsPEV; }
};

template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<3,T,Real>;
    using vector_type = typename base_type::vector_type;
    using matrix_type = typename base_type::matrix_type;

    static constexpr matrix_type _matrix{
      vector_type{5.0,  1.1639777949432226, -0.1639777949432225},
      vector_type{-5.727486121839514,  2.0, 0.7274861218395141},
      vector_type{10.163977794943223,  -9.163977794943223, 5.0}
    };

    static constexpr vector_type _weight{
      Real{0.2777777777777778},
      Real{0.4444444444444444},
      Real{0.2777777777777778} 
    };

    static constexpr vector_type _zero{1,1,1};
    
    static constexpr matrix_type _preconditioner{
      vector_type{0.13888888888888888889, 0, 0},
      vector_type{0.30026319498086459244, 0.22222222222222222222,  0},
      vector_type{0.26798833376246945173, 0.48042111196938334790,  0.13888888888888888889}
    };
    
    static constexpr Real _maxAbsMEV{5.1};
    static constexpr Real _maxAbsPEV{0.22222};

  public:

    static constexpr auto matrix(const size_t i, const size_t j) { return _matrix[j][i]; }
    static constexpr auto weight(const size_t i) { return _weight[i];  }
    static constexpr auto zero(const size_t i) { return _zero[i]; }
    static constexpr auto preconditioner(const size_t i, const size_t j) { return _preconditioner[j][i]; }
    static constexpr auto maxAbsMEV() { return _maxAbsMEV; }
    static constexpr auto maxAbsPEV() { return _maxAbsPEV; }
};


template<size_t ORDER, typename PDE>
struct ODESolverT
{
  using Real           = typename PDE::Real;
  using Vector         = typename PDE::Vector;
  using Expansion      = ExpansionT<ORDER,Vector,Real>;
  using range_iterator = make_range_iterator<size_t>;


  PDE _pde;
  typename Expansion::storage _x, _rhs;

  auto expansionRange() const 
  {
    return make_range_iteratorT<0,Expansion::size()>{};
  }

  ODESolverT(const PDE &pde) : _pde{pde}
  {
    for (auto k : expansionRange())
    {
      _x  [k].resize(_pde.resolution());
      _rhs[k].resize(_pde.resolution());
    }
  };


  void rhs(const Vector &u0)
  {
    for (auto k : expansionRange())
    {
      const auto scale = Expansion::weight(k);
      _pde.compute_rhs(_x[k], _rhs[k], [scale](const auto x) { return x * scale; });

      assert(_x[k].size() == u0.size());
      for (auto i : range_iterator{0,u0.size()})
      {
        Real r = 0;
        for (auto l : expansionRange())
        {
          r += Expansion::matrix(k,l) * (u0[i] - _x[l][i]);
        }
        _rhs[k][i] += r;
      }
    }

    /* precondition */
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

    const Real omega = Real{0.9}/(Expansion::maxAbsPEV()*(Expansion::maxAbsMEV() + _pde.cfl()));
    printf(std::cerr, "omega= %   PEV= %  MEV= %  cfl= % \n",
        omega,
        Expansion::maxAbsPEV(),
        Expansion::maxAbsMEV(),
        _pde.cfl());

    for (auto k : expansionRange())
      for (auto& x : _rhs[k])
        x *= omega;
  }

  void solve_system(const Vector& u0)
  {
    using std::get;
    constexpr auto niter = 10;
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
        for (auto v : make_zip_iterator(_x[k], _rhs[k]))
        {
          get<0>(v) += get<1>(v);
          error[k] += square(get<1>(v)); ///(get<0>(v) + eps));
          cnt += 1;
        }
        err += std::max(err,std::sqrt(error[k]/cnt));
      }
      printf(std::cerr, "iter= %  err= % \n", iter, err);
    }
  }

  void update()
  {
    static auto du = _pde.state();

    for (auto k : expansionRange())
    {
      _x[k] = _pde.state();
    }
    solve_system(_pde.state());

    for (auto i : range_iterator{0,du.size()})
    {
      Real dy = 0;
      du[i] = 0;
      for (auto k : expansionRange())
      {
        du[i] += Expansion::weight(k) * _x[k][i];
      }
    }
    _pde.update(du);
  }
};

template<typename real_type>
struct PDEDiffusion
{
  using Real   = real_type;
  using Vector = std::vector<Real>;
  using range_iterator = make_range_iterator<size_t>;

  Vector _f;

  Real _time;
  Real _cfl;
  Real _dx;
  Real _diff;
  Real _zeta;

  Real dt() const {return _cfl * _diff/square(_dx);}

  auto time() const { return _time; }
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

template<typename PDE>
void dump2file(std::string fileName, const PDE &pde)
{
  const auto &f = pde.state();
  const auto n = f.size();
  const auto dx = pde.dx();
  std::ofstream fout(fileName);
  fout << "# time= " << pde.time() << std::endl;
  fout << "# n= " << n-2 << std::endl;
  for (size_t i = 1; i < n-1; i++)
  {
    const auto x = (i-1+0.5)*dx/(n-2);
    fout << x << " " << std::setprecision(16) << f[i] << std::endl;
  }
}

int main(int argc, char * argv[])
{
  const size_t ncell = argc > 1 ? atoi(argv[1])+2 : 102;
  printf(std::cerr, "ncell= %\n", ncell);

  const size_t niter = argc > 2 ? atoi(argv[2]) : 1;
  printf(std::cerr, "niter= %\n", niter);
  

  using Real = double;

  constexpr auto ORDER = 1;
  using PDE = PDEDiffusion<Real>;
  using Solver = ODESolverT<ORDER,PDE>;

  PDE pde{ncell};

  pde._dx   = 1;
  pde._diff = 1;
  pde._cfl  = 0.40;  /* stable for cfl <= 0.5 */
  pde._time = 0;

  pde.set_ic();
  dump2file("ic.txt",pde);

  Solver solver(pde);

  for (int iter = 0; iter < niter; iter++)
  {
    solver.update();
    printf(std::cerr, "iter= %\n", iter);
    dump2file("iter"+std::to_string(iter),pde);
  }




  return 0;  

}




