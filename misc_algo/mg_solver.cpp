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


    static constexpr auto matrix(const size_t i, const size_t j);
    static constexpr auto maxAbsMEV();
    static constexpr auto preconditioner(const size_t i, const size_t j);
    static constexpr auto maxAbsPEV();

    static constexpr auto weight(const size_t i);
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j);
};

template<size_t N, typename T, typename Real>
class ExpansionT;

template<typename T, typename Real>
class ExpansionT<1,T,Real> : public ExpansionBaseT<1,T,Real>
{
  protected:
    static constexpr size_t N = 1;
    using base_type   = ExpansionBaseT<1,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    {
      return Real{2};
    }
    static constexpr auto maxAbsMEV() 
    {
      return Real{2};
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    { 
      return Real{0.5};
    }
    static constexpr auto maxAbsPEV() 
    { 
      return Real{0.5};
    }

    static constexpr auto weight(const size_t i) 
    {
      return Real{1};
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(2 == ORDER || 3 == ORDER, " Pronlogation order is not supported");
        if (2 == ORDER)
        {
          constexpr Real matrix[2][N]={{1},{1}};
          return matrix[i][j];
        }
        else if (3 == ORDER)
        {
          constexpr Real matrix[3][N]={{1},{1},{1}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  protected:
    static constexpr size_t N = 2;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{3.0,  0.4641016151377546},{-6.464101615137754,  3.0}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV = 3.5;
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = { {0.250000000000000000000,  0},{0.53867513459481288225, 0.250000000000000000000} };
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV = 0.25;
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] = {0.5,0.5};
      return weight[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(3 == ORDER || 4 == ORDER, " Pronlogation order is not supported");
        if (3 == ORDER)
        {
          constexpr Real matrix[3][N]={{1.170820393249936908923,-0.170820393249936908923},{0.500000000000000000000,0.500000000000000000000},{-0.170820393249936908923,1.170820393249936908923}};
          return matrix[i][j];
        }
        else if (4 == ORDER)
        {
          constexpr Real matrix[4][N]={{1.245765921961681556807,-0.245765921961681556807},{0.794432220549629981178,0.205567779450370018822},{0.205567779450370018822,0.794432220549629981178},{-0.245765921961681556807,1.245765921961681556807}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  protected:
    static constexpr size_t N = 3;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{5.0,  1.1639777949432226, -0.1639777949432225},{-5.727486121839514,  2.0, 0.7274861218395141},{10.163977794943223,  -9.163977794943223, 5.0}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{5.1};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = { {0.13888888888888888889, 0, 0}, {0.30026319498086459244, 0.22222222222222222222,  0}, {0.26798833376246945173, 0.48042111196938334790,  0.13888888888888888889}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.22222};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] = {0.2777777777777778,0.4444444444444444,0.2777777777777778};
      return weight[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(4 == ORDER || 5 == ORDER, " Pronlogation order is not supported");
        if (4 == ORDER)
        {
          constexpr Real matrix[4][N]={{1.173824221557882097734,-0.235926245243015346149,0.062102023685133248415},{0.315779411635934322717,0.807354816671586774721,-0.123134228307521097438},{-0.123134228307521097438,0.807354816671586774721,0.315779411635934322717},{0.062102023685133248415,-0.235926245243015346149,1.173824221557882097734}};
          return matrix[i][j];
        }
        else if (5 == ORDER)
        {
          constexpr Real matrix[5][N]={{1.269238169652725404510,-0.368603188642368014803,0.099365018989642610293},{0.589204776685259874958,0.516751336790516162951,-0.105956113475776037910},{0,1.000000000000000000000,0},{-0.105956113475776037910,0.516751336790516162951,0.589204776685259874958},{0.099365018989642610293,-0.368603188642368014803,1.269238169652725404510}};
          return matrix[i][j];
        }
      }
};
template<typename T, typename Real>
class ExpansionT<4,T,Real> : public ExpansionBaseT<4,T,Real>
{
  protected:
    static constexpr size_t N = 4;
    using base_type = ExpansionBaseT<N,T,Real>;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = {{7.73861278752583057,2.045089650303908785,-0.4370708023957989036,0.0866440235032616747},{-7.201340999706890565,2.261387212474169433,1.448782034533681252,-0.2331339812124165654},{6.343622218624971333,-5.971556459482020118,2.261387212474169433,1.090852762294335798},{-15.56386959855492281,11.89278387805684136,-13.50080272596495125,7.73861278752583057}};
      return matrix[i][j]; 
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{6.8};
      return maxAbsMEV;
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = {{0.0869637112843634643,0,0,0},{0.1881181174998680717,0.1630362887156365357,0,0},{0.1671919219741887732,0.3539530060337439665,0.1630362887156365357,0},{0.1774825722545226118,0.3134451147418683468,0.3526767575162718646,0.0869637112843634643}};
      return preconditioner[i][j];
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.17};
      return maxAbsPEV; 
    }

    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] = {0.17392742256872748,0.32607257743127305,0.3260725774312732,0.17392742256872748};
      return weight[i];
    }
    template<size_t ORDER>
      static constexpr auto prolongate(const size_t i, const size_t j)
      {
        static_assert(5==ORDER || 6==ORDER, " Pronlogation order is not supported");
        if (5==ORDER)
        {
          constexpr Real matrix[5][N]={{1.15665233447960049451,-0.23306848457062838863,0.105895713739470568226,-0.029479563648442674104},{0.226361866637914228271,0.93205208186650974176,-0.210599724002141635027,0.052185775497717664992},{-0.092326598440728820911,0.592326598440728820911,0.592326598440728820911,-0.092326598440728820911},{0.052185775497717664992,-0.210599724002141635027,0.93205208186650974176,0.226361866637914228271},{-0.029479563648442674104,0.105895713739470568226,-0.23306848457062838863,1.15665233447960049451}};
          return matrix[i][j];
        }
        else if (6==ORDER)
        {
          constexpr Real matrix[6][N]={{1.25427669651204691779,-0.38249201436315504170,0.178098950270956183941,-0.049883632419848060027},{0.454139585055516588660,0.71591916408312051758,-0.229700084831735121197,0.059641335693098014958},{-0.059826670151450920280,0.93065512687185849208,0.163036458535962670788,-0.033864915256370242588},{-0.033864915256370242588,0.163036458535962670788,0.93065512687185849208,-0.059826670151450920280},{0.059641335693098014958,-0.229700084831735121197,0.71591916408312051758,0.454139585055516588660},{-0.049883632419848060027,0.178098950270956183941,-0.38249201436315504170,1.25427669651204691779}};
          return matrix[i][j];
        }
      }
};

template<size_t ORDER_MAX, typename PDE>
class ODESolverT
{
  public:
    using Real           = typename PDE::Real;
    using Vector         = typename PDE::Vector;
  private:
    using range_iterator = make_range_iterator<size_t>;

    template<size_t ORDER>
      using Expansion = ExpansionT<ORDER,Vector,Real>;

    PDE _pde;
    Real _time;
    bool _verbose;

    Real _err, _err_pre, _cfl_pre, _cfl;
    Real _atol, _rtol;

    typename Expansion<ORDER_MAX>::storage _x, _rhs, _rhs_pde;

    template<size_t K>  
      auto expansionRange() const 
      {
        return make_range_iteratorT<0,K>{};
      }
    
    static constexpr Real omegaCFL = 0.5;

  public:
    auto time() const {return _time;}

    ODESolverT(const PDE &pde) : _pde{pde}, _time{0}
    {
      for (auto k : expansionRange<ORDER_MAX>())
      {
        _x  [k].resize(_pde.resolution());
        _rhs[k].resize(_pde.resolution());
        _rhs_pde[k].resize(_pde.resolution());
        std::fill(_x[k].begin(), _x[k].end(), 0);
      }
      _err = _err_pre = -1;
      _cfl = _cfl_pre = -1;

      const Real tol = 1.0e-10;
      _atol = tol;
      _rtol = tol;
    };
    PDE& pde() { return _pde; }
    const PDE& pde() const { return _pde; }

    template<size_t ORDER>
      void rhs(const Vector &u0)
      {
        using std::get;
        auto x = _x;

        /* compute RHS */
        const auto h = _pde.dt();
        for (auto k : expansionRange<ORDER>())
        {
          for (auto v : make_zip_iterator(x[k], u0))
            get<0>(v) += get<1>(v);
          _pde.compute_rhs(_rhs_pde[k], x[k]);

          assert(_x[k].size() == u0.size());
          for (auto i : range_iterator{u0.size()})
          {
            Real r = 0;
            for (auto l : expansionRange<ORDER>())
              r += Expansion<ORDER>::matrix(k,l) * _x[l][i];
            _rhs[k][i] = h*_rhs_pde[k][i]  - r;
          }
        }

        /* precondition RHS */
        for (auto i : range_iterator{0,u0.size()})
        {
          std::array<Real,Expansion<ORDER>::size()> tmp;
          for (auto k : expansionRange<ORDER>())
          {
            tmp[k] = 0;
            for (auto l : expansionRange<ORDER>())
              tmp[k] += Expansion<ORDER>::preconditioner(k,l)*_rhs[l][i];
          }
          for (auto k : expansionRange<ORDER>())
            _rhs[k][i] = tmp[k];
        }
      }

    template<size_t ORDER>
      void smoother(const int n_smooth_iter, const Vector &u0)
      {
        const auto n = n_smooth_iter;

        using std::get;

        const Real omega = omegaCFL/(Expansion<ORDER>::maxAbsPEV()*(Expansion<ORDER>::maxAbsMEV() + _pde.AbsEV()));

        static decltype(_x) y0, y1,tmp;
        y0  = _x;

        rhs<ORDER>(u0);
        static auto res = _x;


        auto frac = Real{1}/(1+n)/(2+n)/(3+4*n);
        auto omega0 = frac*6;
        auto omega1 = frac*3*n*(3+n);
        auto omegak = [n,frac](const int k) 
        {
          return frac*3*(1-k+n)*(2+k+n);
        };

        for (auto k : expansionRange<ORDER>())
        {
          for (auto v : make_zip_iterator(_x[k], _rhs[k],res[k]))
          {
            auto&   x = get<0>(v);
            auto& rhs = get<1>(v);
            auto& r = get<2>(v);
            r = x*omega0 + (3*x + 4.0*omega*rhs)*omega1;
            x = x + 2.0*omega*rhs;
          }
          y1[k] = _x[k];
        }



        for (int i = 2; i <= n; i++)
        {
          rhs<ORDER>(u0);
          for (auto k : expansionRange<ORDER>())
          {
            for (auto v : make_zip_iterator(_x[k], y1[k], y0[k], _rhs[k],res[k]))
            {
              auto&   x = get<0>(v);
              auto&  y1 = get<1>(v);
              auto&  y0 = get<2>(v);
              auto& rhs = get<3>(v);
              x = 2*y1 - y0 + 4*omega*rhs;
              auto& r = get<4>(v);
              r += 2*omegak(i)*x;
            }
            y0[k] = y1[k];
            y1[k] = _x[k];
          }
        }
        for (auto k : expansionRange<ORDER>())
          _x[k] = res[k];
      }


    template<size_t ORDER_NEW,size_t ORDER_OLD>
      void solve_system_mg(const size_t n_smooth_iter, const Vector &u0)
      {
        const auto n = u0.size();
        const auto h = _pde.dt();

        Vector y0;
        y0.resize(n);

        using ExpansionNew = Expansion<ORDER_NEW>;
        using ExpansionOld = Expansion<ORDER_OLD>;

        /* interpolate from old to new order */
        for (auto i : range_iterator{n})
        {
          std::array<Real,ExpansionNew::size()> x;
          for (auto k : expansionRange<ExpansionNew::size()>())
          {
            x[k] = 0;
            for (auto l : expansionRange<ExpansionOld::size()>())
            {
              const auto weight = ExpansionOld::template prolongate<ExpansionNew::size()>(k,l);
              x[k] += weight * _x[l][i];
            }
          }
          y0[i] = 0;
          for (auto k : expansionRange<ExpansionNew::size()>())
          {
            _x[k][i] = x[k];
            y0[i] += ExpansionNew::weight(k)*_rhs_pde[k][i];
          }
          y0[i] *= h;
        }

        size_t n_iter = 16;

        const Real atol = _atol;
        const Real rtol = _rtol;

        Real err = 0;
        for (auto iter : range_iterator{n_iter})
        {
          smoother<ORDER_NEW>(n_smooth_iter, u0);
          {
            using std::get;
            auto x = _x;
            rhs<ORDER_NEW>(u0);

            /* compute RHS */
            for (auto k : expansionRange<ExpansionNew::size()>())
            {
              for (auto v : make_zip_iterator(x[k], u0))
                get<0>(v) += get<1>(v);
              _pde.compute_rhs(_rhs_pde[k], x[k]);
            }
          }
          for (auto i : range_iterator{n})
          {
            auto du0 = y0[i];
            decltype(du0) du1 = 0;
            for (auto k : expansionRange<ExpansionNew::size()>())
              du1 += ExpansionNew::weight(k)*_rhs_pde[k][i];
            y0[i] = h*du1;

            const auto aerr = std::abs(du1-du0);
            const auto ym  = std::max(std::abs(u0[i]+du0), std::abs(u0[i]+du1));
            err += square(aerr/(atol + rtol*ym));
          }
          err = std::sqrt(err/n);

          if (_verbose)
          {
            printf(std::cerr, " >>  iter= %  err= % \n", iter, err);
          }
          if (err < 1)
            break;

          if (n_iter-1 == iter && _verbose)
            printf(std::cerr, "   ** iter= %  err= % \n ", iter, err);
        }
      }

    void update(const bool verbose = true)
    {
      _verbose = verbose;
      using std::get;
      const auto& u0 = _pde.state();
      const auto n = u0.size();

      const auto h = _pde.dt();

      auto n_smooth_iter = static_cast<size_t>(1+3*std::sqrt(_pde.cfl()));

      for (auto i : range_iterator{n})
      {
        _x      [0][i] = u0[i];
        _rhs_pde[0][i] = 0;
      }

      Vector du_ctrl, du_solv;
      du_ctrl.resize(n);  
      du_solv.resize(n);

      static_assert(4 == ORDER_MAX, " Order mismatch ");

      solve_system_mg<2,1>(n_smooth_iter, u0);
      solve_system_mg<3,2>(n_smooth_iter, u0);
      for (auto i : range_iterator{n})
      {
        du_ctrl[i] = 0;
        for (auto k : expansionRange<ORDER_MAX-1>())
          du_ctrl[i] += Expansion<ORDER_MAX-1>::weight(k)*_rhs_pde[k][i];
        du_ctrl[i] *= h;
      }

      solve_system_mg<4,3>(n_smooth_iter, u0);
      for (auto i : range_iterator{n})
      {
        du_solv[i] = 0;
        for (auto k : expansionRange<ORDER_MAX>())
          du_solv[i] += Expansion<ORDER_MAX>::weight(k)*_rhs_pde[k][i];
        du_solv[i] *= h;
      }


      Real err = 0;
      for (auto i : range_iterator{n})
      {
        const auto u_solv = u0[i] + du_solv[i];
        const auto u_ctrl = u0[i] + du_ctrl[i]; 

        const auto um = std::max(std::abs(u0[i]), std::abs(u_solv));
        const auto sc1 = _atol + _rtol*um;
        const auto du_err = std::abs(u_solv - u_ctrl);
        err += square(du_err/sc1);
      }
      err = std::sqrt(err/n);
      _err_pre  = _err;
      _err      = err;
      _cfl_pre  = _cfl;
      _cfl      = _pde.get_cfl();

      _pde.update(du_solv);
      _time += _pde.dt();

      if (_verbose)
        printf(std::cerr, "# -- err_pre= %   err= %  cfl_pre= %  cfl= %\n", _err_pre, _err,
            _cfl_pre, _cfl);

      Real cfl_scale = 1;
      if (_err > 0 && _err_pre > 0)
      {
        const auto p = Real{1}/(Expansion<ORDER_MAX>::size());
        cfl_scale = 0.8*std::pow(1/_err,p)*_cfl/_cfl_pre*std::pow(_err_pre/_err,p);
      }
      else if (_err > 1)
      {
        const auto p = Real{1}/(Expansion<ORDER_MAX>::size());
        cfl_scale = 0.8*std::pow(1/_err,p);
      }
    }
};

template<typename real_type>
class PDEBurger
{
  public:
    using Real   = real_type;
    using Vector = std::vector<Real>;

  private:
    using range_iterator = make_range_iterator<size_t>;

    Vector _f;

    Real _cfl;
    Real _dx;
    Real _diff;

    size_t n_rhs_calls;

  public:
    void set_dx(const Real dx) { _dx = dx;}
    void set_diff(const Real diff) { _diff = diff;}
    void set_cfl(const Real cfl) { _cfl = cfl;}
    auto get_cfl() const { return _cfl; }
    Real dt() const {
      return  _cfl * 0.5*square(_dx)/_diff;
    }

    auto cost() const { return n_rhs_calls; }
    Real AbsEV() const
    {
      return dt() * 4.0*_diff/square(_dx);  /* 2.0 * cfl */
    }

    auto dx() const { return _dx; }

    PDEBurger(const size_t n) : _f{Vector(n+2)}, n_rhs_calls{0}
    {
    }
    auto cfl() const { return _cfl; }
    auto resolution() const { return _f.size(); }


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
      void compute_rhs(Vector &res, Vector &x, Func func)
      {
        n_rhs_calls++;
//        const auto c = dt() * _diff/square(_dx);
        const auto c =  _diff/square(_dx);
        const auto n = x.size();
        res[0] = c * (x[n-1] - Real{2.0} * x[0] + x[1]);
        res[0] = func(res[0]);
        for (auto i : range_iterator{1,x.size()-1})
        {
          res[i] = c * (x[i-1] - Real{2.0} * x[i] + x[i+1]);
          res[i] = func(res[i]);
        }
        res[n-1] = c * (x[n-2] - Real{2.0} * x[n-1] + x[0]);
        res[n-1] = func(res[n-1]);
      }
    void compute_rhs(Vector &res, Vector &x)
    {
      compute_rhs(res, x, [](const auto x) { return x; });
    }

    void wedge_ic()
    {
      using std::max;

      /* set  a profile of delta function */

      auto &f = _f;

      const int n = f.size();
      const auto dx = _dx;
      const auto L = dx*(n-2);
      const auto dL = L * 0.1;
      const auto ic = n>>1;

      const auto ampl = Real{10.0};
      const auto slope = ampl/dL;


      const auto fmin = Real{1};
      std::fill(f.begin(), f.end(), fmin);
      const int m = static_cast<int>(dL/dx + 0.5);
      for (int i = -m; i <= m; i++)
      {
        const auto x = L/2 + dx*i;
        f[ic - ic/2 + i] = std::max(ampl - slope*(std::abs(L/2-x)),fmin);
      }
      for (int i = -m*2; i <= m*2; i++)
      {
        const auto x = L/2 + dx*i;
        f[ic + ic/2 + i] = std::max(1.5*ampl - slope*(std::abs(L/2-x)),fmin);
      }
    }
    void burger_ic()
    {
      using std::max;

      /* set  a profile of delta function */

      auto &f = _f;

      const int n = f.size();
      const auto dx = _dx;
      for (int i = 0; i < n; i++)
      {
        const auto x = (i+1)*dx;
        f[i] = 1.5*x*square(1-x);
      }
    }
    void set_ic()
    {
    //  wedge_ic();
        burger_ic();
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
    const auto x = (i-1)*dx;
    fout << x << " " << std::setprecision(16) << f[i] << std::endl;
  }
}

  template<typename Solver>
auto compute_mass(const Solver &solver)
{
  const auto &f = solver.pde().state();
  const auto n = f.size();
  const auto dx = solver.pde().dx();
  typename Solver::Real sum = 0;
  for (int i = 0; i < n; i++)
    sum += f[i]*dx;
  return sum;
}

int main(int argc, char * argv[])
{
  using Real = double;

  const size_t ncell = argc > 1 ? atoi(argv[1]) : 128;
  printf(std::cerr, "ncell= %\n", ncell);

  const Real tend = argc > 2 ? atof(argv[2]) : 0.005;
  printf(std::cerr, "tend= %\n", tend);
  


  constexpr auto ORDER = 4;
  using PDE = PDEBurger<Real>;
  using Solver = ODESolverT<ORDER,PDE>;


  Solver solver(PDE{ncell});

  solver.pde().set_dx(1.0/ncell);
  solver.pde().set_diff(1);
  solver.pde().set_cfl(0.8);//*2*2*2); //*2*2); //*8); //*8); //*64); //*16); //*64); //*64); //*64); //*64); //*64);//*64); //*64);//*64); //*4); //*64/4); //*64); //*64); //*64/4); //*64*4);//*64); //*64); //*64); //*4*4*4);  /* stable for cfl <= 0.5 */

  const auto dt = solver.pde().dt();

  solver.pde().set_ic();
  const auto mass0 = compute_mass(solver);
  dump2file(solver, "ic.txt");
  {
    size_t kstep = 0;
    bool keep_looping = true;
    while (keep_looping)
    {
      auto verbose_step = (kstep%1) == 0;
      auto verbose_iter = (kstep%1) == 0;
      const auto dt = solver.pde().dt();
      bool break_loop = true;

      if (solver.time() + dt >= tend)
      {
        const auto dt_new = tend - solver.time();
        assert(dt_new >= 0);
        solver.pde().set_cfl(solver.pde().get_cfl() * dt_new/dt);
        keep_looping = false;
        verbose_step = verbose_iter = true;
      }

      if (verbose_step)
      {
        const auto mass = compute_mass(solver);
        printf(std::cerr, "step= % : time= % dt= % (cfl= %) ORDER= % cost= % -- mass_err= %  Tend= % \n", 
            kstep, solver.time(), solver.pde().dt(), solver.pde().cfl(), ORDER,
            solver.pde().cost(), (mass-mass0)/mass0, tend);
      }
//      if (kstep > 45)
//        solver.pde().set_cfl(1000);
      solver.update(verbose_iter);
      kstep++;
    }
#if 0
    const auto mass = compute_mass(solver);
    printf(std::cerr, "step= % : time= % dt= % (cfl= %) ORDER= % cost= % -- mass_err= %  Tend= % \n", 
        kstep, solver.time(), solver.pde().dt(), solver.pde().cfl(), ORDER,
        solver.pde().cost(), (mass-mass0)/mass0, tend);
#endif
  }
  printf(std::cerr, " Writing output ... \n");
//  dump2file(solver);
  dump2file(solver, "burger");
  printf(std::cerr, "cost= %\n", solver.pde().cost());



  return 0;  

}




