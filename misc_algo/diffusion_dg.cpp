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

#if 1
#define PIF_
#else
#define DGF_
#endif

#ifdef PIFT_
template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  /* PIF */
  protected:
    using base_type  = ExpansionBaseT<3,T,Real>;
    static constexpr size_t N = 3;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {11.464101615137755,  0.6188021535170061, -0.0829037686547607},
        {-8.618802153517006,  2., 0.6188021535170061},
        {16.08290376865476, -8.618802153517006, 4.535898384862246}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.22222222222222274,
        0.5555555555555555,
        0.22222222222222274
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.075026719286759500829, 0,0},
        {0.25544867840851755224,0.27777777777777777778, 0},
        {0.21936428658416594196,0.56645291237259066003, 0.147195502935462721393}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{10.8};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.28};
      return maxAbsPEV; 
    }
};
#elif defined PIF_  /* not PIFt_ */
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
      return Real{0.5};
    }
};
template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  /* PIF */
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
  /* PIF */
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
template<typename T, typename Real>
class ExpansionT<5,T,Real> : public ExpansionBaseT<5,T,Real>
{
  /* PIF */
  protected:
    static constexpr size_t N = 5;
    using base_type  = ExpansionBaseT<N,T,Real>;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {11.18330013267037774,3.13131216201181084,-0.758731795980807391,0.239101223353686050,-0.0543147605653389237},
        {-9.44759960151614989,2.81669986732962226,2.21788633227481812,-0.557122620293797301,0.118357949604667387},
        {6.42011650355933662,-6.22012045466975167,2.00000000000000000,1.86599528883177949,-0.315991337721364445},
        {-8.01592078481097030,6.19052235495304182,-7.39311627638238017,2.81669986732962226,1.55003676630984697},
        {22.4209150259060944,-16.1933902399992350,15.4154432215698509,-19.0856011786573598,11.18330013267037774}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.11846344252809503,
        0.2393143352496833,
        0.28444444444444444,
        0.23931433524968349,
        0.1184634425280951
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.0592317212640472719,0,0,0,0},
        {0.128151005670045283,0.119657167624841617,0,0,0},
        {0.1137762880042246025,0.260004651680641519,0.142222222222222222,0,0},
        {0.121232436926864147,0.228996054578999877,0.309036559064086645,0.119657167624841617,0},
        {0.1168753295602285452,0.244908128910495419,0.273190043625801489,0.258884699608759272,0.0592317212640472719}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{8.6};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.15};
      return maxAbsPEV; 
    }
};
#elif defined DGF_  /* not PIF_ */
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
      return Real{1};
    }
    static constexpr auto maxAbsMEV() 
    {
      return Real{1};
    }
    static constexpr auto maxAbsPEV() 
    { 
      return Real{1};
    }
};
template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  /* DGF */
  protected:
    using base_type = ExpansionBaseT<2,T,Real>;
    static constexpr size_t N = 2;

  public:
    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {2.0,  0.7320508075688785},
        {-2.732050807568876,  1.9999999999999958}
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
        {0.3333333333333329,  0},
        {0.4553418012614798,  0.33333333333333365}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV = 2.5;
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV = 0.33;
      return maxAbsPEV; 
    }
};
template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  /* DGF */
  protected:
    using base_type  = ExpansionBaseT<3,T,Real>;
    static constexpr size_t N = 3;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {4.000000000000003, 1.6147844564602516, -0.2909944487358082},
        {-3.509240285287669,  1.0000000000000007, 1.009240285287668},
        {2.290994448735803, -5.6147844564602405, 3.9999999999999774}
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
        {0.16111111111111154, 0,0},
        {0.2724854172030868,  0.2777777777777769, 0},
        {0.29021055598469203, 0.43597666752493847,  0.16111111111111140}
      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{4.1};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.28};
      return maxAbsPEV; 
    }
};
template<typename T, typename Real>
class ExpansionT<4,T,Real> : public ExpansionBaseT<4,T,Real>
{
  /* DGF */
  protected:
    static constexpr size_t N = 4;
    using base_type  = ExpansionBaseT<N,T,Real>;

  public:

    static constexpr auto matrix(const size_t i, const size_t j) 
    { 
      constexpr Real matrix[N][N] = 
      {
        {6.738612787525834, 2.5779942727073806, -0.6995574654955786,  0.1612563383245323},
        {-5.324832600342981,  1.2613872124741736, 1.941340462561435,  -0.3731446166705182},
        {2.5339048478804855,  -3.9413404625614348,  1.2613872124741656, 1.3751045941403923},
        {-2.161256338324527,  4.750469319393827,  -9.982795494470142, 6.738612787525786}
      };
      return matrix[j][i]; 
    }
    static constexpr auto weight(const size_t i)  
    { 
      constexpr Real weight[] =
      {
        0.17392742256872748, 
        0.32607257743127305, 
        0.3260725774312732, 
        0.17392742256872748
      };
      return weight[i];
    }
    static constexpr auto preconditioner(const size_t i, const size_t j) 
    {
      constexpr Real preconditioner[N][N] = 
      {
        {0.0950400941860567,  0,0,0},
        {0.17720653136163217, 0.19067419152822931,  0,0},
        {0.17810350811242565, 0.32631510322115087,  0.19067419152822807,  0},
        {0.1694061893528294,  0.3339017452341196, 0.33222012702401943,  0.09504009418605698}

      };
      return preconditioner[j][i];
    }
    static constexpr auto maxAbsMEV() 
    {
      constexpr Real maxAbsMEV{5.8};
      return maxAbsMEV;
    }
    static constexpr auto maxAbsPEV() 
    { 
      constexpr Real maxAbsPEV{0.20};
      return maxAbsPEV; 
    }
};
#endif


template<size_t ORDER, typename PDE>
class ODESolverT
{
  public:
    using Real           = typename PDE::Real;
    using Vector         = typename PDE::Vector;
  private:
    using Expansion      = ExpansionT<ORDER,Vector,Real>;
    using range_iterator = make_range_iterator<size_t>;


    PDE _pde;
    Real _time;
    bool _verbose;
    typename Expansion::storage _x, _rhs;

    static constexpr Real omegaCFL = 0.9;

    std::array<bool,Expansion::size()> updateExpansion;

  public:
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
      auto x = _x;
      auto& rhs = _rhs;

      /* compute RHS */
      for (auto k : expansionRange())
        if (updateExpansion[k])
        {
          for (auto v : make_zip_iterator(x[k], u0))
            get<0>(v) += get<1>(v);
          _pde.apply_bc(x[k]);
          _pde.compute_rhs(rhs[k], x[k]);

          assert(_x[k].size() == u0.size());
          for (auto l : expansionRange())
            for (auto v : make_zip_iterator(rhs[k], _x[l]))
            {
              get<0>(v) += - Expansion::matrix(k,l) * get<1>(v);
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

    void iterate(const Vector &u0, int n)
    {
      using std::get;

      const Real omega = omegaCFL/(Expansion::maxAbsPEV()*(Expansion::maxAbsMEV() + _pde.AbsEV()));

      static decltype(_x) y0, y1,tmp;
      y0  = _x;

      rhs(u0);
      static auto res = _x;

#if 1
#define WP   /* prefered for high resolution & large time step, use scaleing 1x for nstate with cfl */
#elif 1
#define OPT
#endif

#ifdef WP
      auto scale = Real{1}/(2*n+1);
#elif defined OPT
      auto frac = Real{1}/(1+n)/(2+n)/(3+4*n);
      auto omega0 = frac*6;
      auto omega1 = frac*3*n*(3+n);
      auto omegak = [n,frac](const int k) 
      {
        return frac*3*(1-k+n)*(2+k+n);
      };
#endif

      for (auto k : expansionRange())
        for (auto v : make_zip_iterator(_x[k], _rhs[k],res[k]))
        {
          auto&   x = get<0>(v);
          auto& rhs = get<1>(v);
#ifdef OPT
          auto& r = get<2>(v);
          r = x*omega0 + (3*x + 4.0*omega*rhs)*omega1;
#elif defined WP
          auto& r = get<2>(v);
          r = (3*x + 4.0*omega*rhs)*scale;
#else
          x = x + 2.0*omega*rhs;
#endif
        }


      y1 = _x;

      for (int i = 2; i <= n; i++)
      {
        rhs(u0);
        for (auto k : expansionRange())
          for (auto v : make_zip_iterator(_x[k], y1[k], y0[k], _rhs[k],res[k]))
          {
            auto&   x = get<0>(v);
            auto&  y1 = get<1>(v);
            auto&  y0 = get<2>(v);
            auto& rhs = get<3>(v);
#ifdef OPT
            x = 2*y1 - y0 + 4*omega*rhs;
            auto& r = get<4>(v);
            r += 2*omegak(i)*x;
#elif defined WP
            x = 2*y1 - y0 + 4*omega*rhs;
            auto& r = get<4>(v);
            r += 2*scale*x;
#else
            const auto a = (Real{2}*i-1)/i;
            const auto b = Real{-1}*(i-1)/i;
            x = a*y1 + b*y0 + 2*omega*a*rhs;
#endif
          }
        y0 = y1;
        y1 = _x;
      }
#ifdef OPT
      _x = res;
#elif defined WP
      _x = res;
#endif

    }


    void iterate(const Vector &u0, bool verbose)
    {
      const int nstage = static_cast<int>(1+std::sqrt(_pde.cfl()));  /* stiffff */
      iterate(u0, nstage);
      if (verbose)
      {
        printf(std::cerr, " nstage= % \n", nstage);
      }
    }

    void solve_system(const Vector& u0)
    {
      using std::get;
      size_t  niter = 5; //8*2*2; // * 32; //*2; //16 ;//1; //32; //50;
      niter = 31;
      std::array<Real,Expansion::size()> error;
      constexpr Real tol = 1.0e-5;
      bool verbose = _verbose;
      for (auto k : expansionRange())
        updateExpansion[k] = true;

      for (auto iter : range_iterator{0,niter})
      {
        for (auto k : expansionRange())
        {
          error[k] = 0;
        }

        iterate(u0, verbose);
        verbose = false;

        Real err = 0;
        int cnt = 0;
        constexpr auto eps = Real{1.0e-7};
        for (auto k : expansionRange())
        {
          for (auto i : range_iterator{1,u0.size()-1})
          {
            auto err = square(_rhs[k][i]);
    //        err *= 1.0/square(_x[k][i] + eps);
            error[k] += err;
         //   error[k] = std::max(error[k], std::abs(_rhs[k][i]));
            cnt += 1;
          }
          error[k] = std::sqrt(error[k]);
          error[k] *= std::sqrt(Real{1}/cnt);
          err = std::max(err,error[k]);
        }
        if (_verbose)
        {
          printf(std::cerr, " >>  iter= %  ", iter);
          for (auto k : expansionRange())
          {
            updateExpansion[k] = error[k] > tol;
            printf(std::cerr, "e[%]= % ", k, error[k]);
          }
          printf(std::cerr, "\n");
        }
        if (err <= tol)
        {
          if (_verbose)
          {
            printf(std::cerr, "   >> iter= %  err= % \n", iter, err);
          }
          break;
        }
        if (iter == niter - 1 && _verbose)
          printf(std::cerr, "   ** iter= %  err= % \n", iter, err);
      }
    }

    void update(const bool verbose = true)
    {
      _verbose = verbose;
      using std::get;
      static auto du = _pde.state();

      for (auto k : expansionRange())
        for (auto &x : _x[k])
          x = 0;

      solve_system(_pde.state());

      for (auto k : expansionRange())
      {
        for (auto v : make_zip_iterator(_x[k], _pde.state()))
          get<0>(v) += get<1>(v);
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
class PDEDiffusion
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
    Real dt() const {return _cfl * 0.5*square(_dx)/_diff;}

    auto cost() const { return n_rhs_calls; }
    Real AbsEV() const
    {
      return dt() * 4.0*_diff/square(_dx);  /* 2.0 * cfl */
    }

    auto dx() const { return _dx; }

    PDEDiffusion(const size_t n) : _f{Vector(n)}, n_rhs_calls{0}
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
      periodic_bc(f);
      //free_bc(f);
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
        n_rhs_calls++;
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

int main(int argc, char * argv[])
{
  using Real = double;

  const size_t ncell = argc > 1 ? atoi(argv[1]) : 128;
  printf(std::cerr, "ncell= %\n", ncell);

  const Real tau = argc > 2 ? atof(argv[2]) : 0.005;
  printf(std::cerr, "tau= %\n", tau);
  


  constexpr auto ORDER = 3;
  using PDE = PDEDiffusion<Real>;
  using Solver = ODESolverT<ORDER,PDE>;


  Solver solver(PDE{ncell+2});

  solver.pde().set_dx(1.0/ncell);
  solver.pde().set_diff(1);
  solver.pde().set_cfl(0.8*64*64); //*64*4);//*64); //*64); //*64); //*4*4*4);  /* stable for cfl <= 0.5 */

  const auto dt = solver.pde().dt();
  const size_t nstep = 1 + std::max(size_t{0}, static_cast<size_t>(tau/dt));

  printf(std::cerr, "dt= %  tau= %   nstep= %\n",
      dt, tau, nstep);

  solver.pde().set_ic();
  dump2file(solver, "ic.txt");
  for (size_t step = 1; step <= nstep; step++)
  {
    auto verbose_step = (step-1)%1 == 0;
    auto verbose_iter = (step-1)%1  == 0;
    if (step == nstep || step == 1)
      verbose_step = verbose_iter = true;
    solver.update(verbose_iter);
    if (verbose_step)
      printf(std::cerr, "step= % : time= %  cost= %\n", step, solver.time(),
          solver.pde().cost());
  }
  printf(std::cerr, " Writing output ... \n");
  dump2file(solver);
  printf(std::cerr, "cost= %\n", solver.pde().cost());



  return 0;  

}




