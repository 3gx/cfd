#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "printf.h"
#include "zip_iterator.h"
#include "common.h"

template<size_t N, typename T, typename Real>
class ExpansionBaseT
{
  protected:
    using vector_type = std::array<Real,N>;
    using matrix_type = std::array<vector_type,N>;
    std::array<T,N> data;
  public:
    using value_type = T;
    static constexpr auto size() { return N; }
    T& operator[](const size_t i) {return data[i];}
    const T& operator[](const size_t i) const {return data[i];}
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

    static constexpr matrix_type _matrix{vector_type{1}};
    static constexpr matrix_type _matrix_inv{vector_type{1}};
    static constexpr vector_type  _weight{1};
    static constexpr vector_type  _zero{1};

  public:
    static constexpr auto matrix(const size_t i, const size_t j) { return _matrix[j][i]; }
    static constexpr auto matrix_inv(const size_t i, const size_t j) { return _matrix_inv[j][i]; }
    static constexpr auto weight(const size_t i) {  return _weight[i];  }
    static constexpr auto zero(const size_t i) { return _zero[i]; }
};

template<typename T, typename Real>
class ExpansionT<2,T,Real> : public ExpansionBaseT<2,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<2,T,Real>;
    using vector_type = typename base_type::vector_type;
    using matrix_type = typename base_type::matrix_type;

    static constexpr matrix_type _matrix{
      vector_type{Real{1},Real{0.36602540378443865}},
        vector_type{Real{-1.3660254037844386},Real{1}} };
    static constexpr matrix_type _matrix_inv{
      vector_type{0.6666666666666666, -0.24401693585629244},
        vector_type{0.9106836025229591, 0.6666666666666666} };
    static constexpr vector_type  _weight{Real{0.5},Real{0.5}};
    static constexpr vector_type  _zero{Real{1.3660254037844386}, Real{-0.36602540378443865}};

  public:
    static constexpr auto matrix(const size_t i, const size_t j) { return _matrix[j][i]; }
    static constexpr auto matrix_inv(const size_t i, const size_t j) { return _matrix_inv[j][i]; }
    static constexpr auto weight(const size_t i)  { return _weight[i];  }
    static constexpr auto zero(const size_t i) { return _zero[i]; }
};

template<typename T, typename Real>
class ExpansionT<3,T,Real> : public ExpansionBaseT<3,T,Real>
{
  protected:
    using base_type   = ExpansionBaseT<3,T,Real>;
    using vector_type = typename base_type::vector_type;
    using matrix_type = typename base_type::matrix_type;

    static constexpr matrix_type _matrix{
      vector_type{ 1.1111111111111112,  0.4485512379056266, -0.08083179131550157},
        vector_type{-1.5596623490167376,  0.4444444444444444,  0.4485512379056266},
        vector_type{ 0.6363873468710571, -1.5596623490167376,  1.1111111111111112} };

    static constexpr matrix_type _matrix_inv{
      vector_type{0.58, -0.18094750193111253, 0.11524199845510998},
        vector_type{0.9809475019311126, 0.625, -0.18094750193111253},
        vector_type{1.04475800154489, 0.9809475019311126, 0.58} };

    static constexpr vector_type  _weight{
      Real{0.2777777777777778},
        Real{0.4444444444444444},
        Real{0.2777777777777778} };

    static constexpr vector_type  _zero{
      Real{1.4788305577012362},
        Real{-0.6666666666666666},
        Real{0.18783610896543051} };

  public:

    static constexpr auto matrix(const size_t i, const size_t j) { return _matrix[j][i]; }
    static constexpr auto matrix_inv(const size_t i, const size_t j) { return _matrix_inv[j][i]; }
    static constexpr auto weight(const size_t i) { return _weight[i];  }
    static constexpr auto zero(const size_t i) { return _zero[i]; }
};


template<typename Real>
struct ParamT
{
  using value_type = Real;
  Real time;
  Real cfl;
  Real dx;
  Real diff;

  Real dt() const {return cfl * diff/square(dx);}
};

struct PeriodicBC
{
  template<typename Vector>
    static void apply(Vector &f) 
    {
      const auto n = f.size();
      f[0  ] = f[n-2];
      f[n-1] = f[1  ];
    }
};

struct FreeBC 
{
  template<typename Vector>
    static void apply(Vector &f) 
    {
      const auto n = f.size();
      f[0  ] = f[  1];
      f[n-1] = f[n-2];
    }
};

template<typename Param, typename Vector, typename Func>
static void compute_df(const Param &params, const Vector &f, Vector &df, const Func &func)
{
  using Real = typename Vector::value_type;
  const auto c = params.diff/square(params.dx);
  const int n = f.size();
  for (int i = 1 ; i < n-1; i++)
  {
    df[i] = c * (f[i+1] - Real{2.0} * f[i] + f[i-1]);
    df[i] = func(df[i]);
  }
}

  template<typename Param, typename Vector>
static void compute_df(const Param &params, const Vector &f, Vector &df)
{
  compute_df(params, f, df, [](const auto x) { return x; });
}

template<typename Expansion, typename Param, typename PDE>
static void compute_rhs_preconditioned(
    Expansion &rhs,
    const Expansion &x, 
    const typename Expansion::value_type &x0, 
    const Param &param,
    const PDE &pde)
{
  constexpr auto expansionOrder = Expansion::size();   // expansion order

  const auto arraySize = x[0].size();

  using ExpansionType = typename Expansion::value_type;
  using Real = typename ExpansionType::value_type;

  static ExpansionType diff_x;

  for (int k = 0; k < expansionOrder; k++)
  {
    const auto weight = Expansion::weight(k);
    const auto scale = weight * param.dt();
    pde(param, x[k], diff_x, [scale](const auto &val) { return val*scale;} );
    for (int i = 0; i < arraySize; i++)
    {
      Real r = 0;
      for (int l = 0; l < expansionOrder; l++)
        r += Expansion::matrix(k,l)*x[l][i];
      rhs[k][i] = param.zeta*(-r + Expansion::zero(k)*x0[i] + diff_x[i]);
    }
  }
}

template<typename BC, typename Expansion, typename Param, typename PDE>
static void solve_system(
    Expansion &x, 
    const typename Expansion::value_type &x0, 
    const Param &param,
    const PDE &pde)
{
  using ExpansionType = typename Expansion::value_type;
  using Real = typename ExpansionType::value_type;
  using std::get;

  constexpr auto expansionOrder = Expansion::size();   // expansion order

  Expansion rhs;
  const int niter = 10;
  for (int iter = 0; iter < niter; iter++)
  {
    for (int k = 0; k < expansionOrder; k++)
      BC::apply(x[k]);

    iterate(rhs, x, x0, param, pde);
    for (int k = 0; k < expansionOrder; k++)
    {
      Real err = 0;
      constexpr Real eps = 1.0e-7;
      for (auto v : make_zip_iterator(x[k], rhs[k]))
      {
        get<0>(v) += get<1>(v);
        err += square(get<1>(v)/(get<0>(v) + eps));
      }
    }
  }
}


template<size_t ORDER, typename BC, typename Vector, typename Param, typename PDE>
static void update_step(
    Vector &y,
    const Param &param,
    const PDE &pde)
{
  using Real = typename Vector::value_type;
  using Expansion = ExpansionT<ORDER, Vector, Real>;

  Expansion x;
  constexpr int expansionOrder = Expansion::size();
  for (auto k = 0*expansionOrder; k < expansionOrder; k++)
  {
    x[k] = y;
  }

  solve_system<BC>(x, y, param, pde);

  const auto arraySize = y.size();
  const auto dt = param.dt();
  for (auto i = 0*arraySize; i < arraySize; i++)
  {
    Real r = 0;
    for (auto k = 0*expansionOrder; k < expansionOrder; k++)
      r += Expansion::weight(k)*x[k][i];
    y[i] += param.dt*r;
  }
}


template<typename BC,typename Param, typename Vector>
static void compute_update(const Param &params, Vector &f)
{
  using Real = typename Vector::value_type;
  const auto n = f.size();

  BC::apply(f);
  

  static Vector df(n);

  compute_df(params,f,df);

  const auto dt = params.dt();

  for (int i = 1; i < n-1; i++)
  {
    f[i] += dt * df[i];
  }

  BC::apply(f);
}

template<typename BC, typename Param, typename Vector>
static void compute_update_dg1(const Param &params, Vector &f)
{
  using Real = typename Vector::value_type;
  const auto n = f.size();
  const auto dt = params.dt();

  const Vector f0 = f;
  const int niter = 10;
  for (int iter = 0; iter < niter; iter++)
  {
    static Vector df(n);
    BC::apply(f);
    compute_df(params,f,df);
    for (int i = 1; i < n-1; i++)
    {
      f[i] = f0[i] + dt * df[i];
    }
  }

  BC::apply(f);

}

template<typename Param, typename Vector>
void set_ic(const Param &params, Vector &f)
{
  using Real = typename Vector::value_type;
  using std::max;

  /* set  a profile of delta function */

  const int n = f.size();
  const auto dx = params.dx;
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

template<typename Vector>
void dump2file(const Vector &f, const std::string fileName)
{
  const int n = f.size();
  std::ofstream fout(fileName);
  for (int i = 1; i < n-1; i++)
    fout << std::setprecision(16) << f[i] << std::endl;
}

int main(int argc, char * argv[])
{
  const int ncell = argc > 1 ? atoi(argv[1]) : 100;
  printf(std::cerr, "ncell= %\n", ncell);

  const int niter = argc > 2 ? atoi(argv[2]) : 10;
  printf(std::cerr, "niter= %\n", niter);
  

  using Real = double;

  using Param = ParamT<Real>;
  using Vector = std::vector<Real>;

  auto params = Param{};
  params.dx   = 1;
  params.diff = 1;
  params.cfl  = 0.40;  /* stable for cfl <= 0.5 */
  params.time = 0;
//  params.cfl  = 0.252;

  Vector f(ncell+2);

  set_ic(params,f);
  dump2file(f,"ic.txt");

  for (int iter = 0; iter < niter; iter++)
  {
    printf(std::cerr, "iter= %\n", iter);
#if 1
    compute_update<FreeBC>(params,f);
//    compute_update<PeriodicBC>(params,f);
#elif 1
    compute_update_dg1<FreeBC>(FreeBC,params,f);
#endif
    dump2file(f,"iter"+std::to_string(iter));
  }




  return 0;  

}




