#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "printf.h"
#include "zip_iterator.h"
#include "common.h"

template<size_t N, typename Real>
struct DGMatrixBase
{
  using value_type = Real;
  using array = std::array<Real,N>;
  using matrix = std::array<array,N>;
  static constexpr auto size() { return N; }
};


template<size_t N, typename Real> struct DGMatrix;

template<typename Real>
struct DGMatrix<1,Real> : public DGMatrixBase<1,Real>
{
  using base   = DGMatrixBase<1,Real>;
  using array  = typename base::array;
  using matrix = typename base::matrix;
  static constexpr matrix _data{array{1}};
  static constexpr auto getm(const size_t i, const size_t j) { return _data[j][i]; }

  static constexpr matrix _matrix_inv{array{1}};
  static constexpr auto getinv(const size_t i, const size_t j) { return _matrix_inv[j][i]; }

  static constexpr array  _weight{1};
  static constexpr auto getw(const size_t i)                 { return _weight[i];  }
};

template<typename Real>
struct DGMatrix<2,Real> : public DGMatrixBase<2,Real>
{
  using base   = DGMatrixBase<2,Real>;
  using array  = typename base::array;
  using matrix = typename base::matrix;
  static constexpr matrix _data{
    array{Real{1},Real{0.36602540378443865}},
    array{Real{-1.3660254037844386},Real{1}}
  };
  static constexpr auto getm(const size_t i, const size_t j) { return _data[j][i]; }

  static constexpr matrix _matrix_inv{
     array{0.6666666666666666, -0.24401693585629244},
     array{0.9106836025229591, 0.6666666666666666}
  };
  static constexpr auto getinv(const size_t i, const size_t j) { return _matrix_inv[j][i]; }

  static constexpr array  _weight{Real{0.5},Real{0.5}};
  static constexpr auto getw(const size_t i)                 { return _weight[i];  }

  static constexpr array  _zerovec{Real{1.3660254037844386}, Real{-0.36602540378443865}};
  static constexpr auto getzero(const size_t i) { return _zerovec[i]; }
};

template<typename Real>
struct DGMatrix<3,Real> : public DGMatrixBase<3,Real>
{
  using base   = DGMatrixBase<2,Real>;
  using array  = typename base::array;
  using matrix = typename base::matrix;
  static constexpr matrix _data{
    array{ 1.1111111111111112,  0.4485512379056266, -0.08083179131550157},
    array{-1.5596623490167376,  0.4444444444444444,  0.4485512379056266},
    array{ 0.6363873468710571, -1.5596623490167376,  1.1111111111111112}
  };
  static constexpr auto getm(const size_t i, const size_t j) { return _data[j][i]; }

  static constexpr matrix _matrix_inv{
     array{0.58, -0.18094750193111253, 0.11524199845510998},
     array{0.9809475019311126, 0.625, -0.18094750193111253},
     array{1.04475800154489, 0.9809475019311126, 0.58}
  };
  static constexpr auto getinv(const size_t i, const size_t j) { return _matrix_inv[j][i]; }

  static constexpr array  _weight{
    Real{0.2777777777777778},
    Real{0.4444444444444444},
    Real{0.2777777777777778}
  };
  static constexpr auto getw(const size_t i)                 { return _weight[i];  }


  static constexpr array  _zerovec{
    Real{1.4788305577012362},
    Real{-0.6666666666666666},
    Real{0.18783610896543051}
  };
  static constexpr auto getzero(const size_t i) { return _zerovec[i]; }
};

#if 0
template<typename Real>
struct DGMatrix<3,Real> : public DGMatrixBase<3,REAL>
{
  static consttexpr matrix _data{
    array{Real{1},Real{0.36602540378443865}},
    array{Real{-1.3660254037844386},Real{1}}
  };
  static constexpr auto get(const size_t i, const size_t j) { return _data[j][i]; }
};
#endif

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

template<typename Param, typename Vector>
static void compute_df(const Param &params, const Vector &f, Vector &df)
{
  using Real = typename Vector::value_type;
  const auto c = params.diff/square(params.dx);
  const int n = f.size();
  for (int i = 1 ; i < n-1; i++)
  {
    df[i] = c * (f[i+1] - Real{2.0} * f[i] + f[i-1]);
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




