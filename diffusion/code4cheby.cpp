# include <iostream>
#include <vector>
#include <fstream>
#include <string>

template<typename T>
static T square(const T x) { return x*x; }

static void fprintf(std::ostream &os, const char *s)
{
  while (*s) {
    if (*s == '%') {
      if (*(s + 1) == '%') {
        ++s;
      }
      else {
        throw "invalid format string: missing arguments";
      }
    }
    os << *s++;
  }
}

template<typename T, typename... Args>
static void fprintf(std::ostream &os, const char *s, T& value, Args... args)
{
  while (*s) {
    if (*s == '%') {
      if (*(s + 1) == '%') {
        ++s;
      }
      else {
        os << value;
        fprintf(os, s + 1, args...); // call even when *s == 0 to detect extra arguments
        return;
      }
    }
    os << *s++;
  }
  throw "extra arguments provided to printf";
}

template<typename real_t>
struct Params
{
  using value_type = real_t;
  real_t diff;
  real_t dx;
  real_t cfl;

  real_t dt() const {return cfl * diff/square(dx);}
};

struct PeriodicBC
{
  template<typename vector_t>
    static void apply(vector_t &f) 
    {
      const auto n = f.size();
      f[0  ] = f[n-2];
      f[n-1] = f[1  ];
    }
};

struct FreeBC 
{
  template<typename vector_t>
    static void apply(vector_t &f) 
    {
      const auto n = f.size();
      f[0  ] = f[  1];
      f[n-1] = f[n-2];
    }
};

template<typename param_t, typename vector_t>
static void compute_df(const param_t &params, const vector_t &f, vector_t &df)
{
  using real_t = typename vector_t::value_type;
  const auto c = params.diff/square(params.dx);
  const int n = f.size();
  for (int i = 1 ; i < n-1; i++)
  {
    df[i] = c * (f[i+1] - real_t{2.0} * f[i] + f[i-1]);
  }
}

template<typename BC,typename param_t, typename vector_t>
static void compute_update(const param_t &params, vector_t &f)
{
  using real_t = typename vector_t::value_type;
  const auto n = f.size();

  BC::apply(f);
  

  static vector_t df(n);

  compute_df(params,f,df);

  const auto dt = params.dt();

  for (int i = 1; i < n-1; i++)
  {
    f[i] += dt * df[i];
  }

  BC::apply(f);
}

template<typename BC, typename param_t, typename vector_t>
static void compute_update_cheby(const param_t &params, vector_t &f, int niter = 10)
{
  using real_t = typename vector_t::value_type;
  const auto n = f.size();
  const auto dt = params.dt();


  const auto cfl = params.cfl; //diff*dt/square(params.dx);
  const auto eigen_min = real_t{0.5};
  const auto eigen_max = cfl*real_t{4.0};

  const auto c_coeff = (eigen_max - eigen_min) * real_t{0.5};
  const auto d_coeff = (eigen_max + eigen_min) * real_t{0.5};

  const auto am = -cfl;
  const auto ac = real_t{1}+real_t{2}*cfl;
  const auto ap = -cfl;


  BC::apply(f);

  const auto f0 = f;
  auto r = f;
  for (int i = 1; i < n-1; i++)
    r[i] = f0[i] - (am*f[i-1]+ac*f[i]+ap*f[i+1]);

  BC::apply(r);

  auto p = r;
  auto alpha = real_t{0.0};
  for (int iter = 0; iter < niter; iter++)
  {
    for (int i = 1; i < n-1; i++)
    {
      if (iter == 0)
      {
        p[i] = r[i];
        alpha = real_t{2.0}/d_coeff;
      }
      else
      {
        const auto beta = square(c_coeff*alpha*real_t{0.5});
        alpha = real_t{1.0}/(d_coeff-beta);
        p[i] = r[i] + beta*p[i];
      }
    }

    BC::apply(p);
    for (int i = 1; i < n-1; i++)
    {
      const auto qi = am*p[i-1] + ac*p[i] + ap*p[i+1];
      f[i] += alpha * p[i];
      r[i] -= alpha * qi;
    }
    BC::apply(f);
    BC::apply(r);
  }


  BC::apply(f);

}

template<typename param_t, typename vector_t>
void set_ic(const param_t &params, vector_t &f)
{
  using real_t = typename vector_t::value_type;
  using std::max;

  /* set  a profile of delta function */

  const int n = f.size();
  const auto dx = params.dx;
  const auto L = dx*(n-2);
  const auto ic = n>>1;

  const auto ampl = real_t{1.0};


  const int di = 10;
  std::fill(f.begin(), f.end(), 0);
  for (int i = -di; i <= di; i++)
  {
    const auto fi = max(ampl*(di - abs(i)),real_t{0});
    f[ic + i] = fi;
  }
  std::fill(f.begin(), f.end(), 0);
  for (int i = -di; i <= di; i++)
  {
    const auto fi = max(ampl*(di - abs(i)),real_t{0});
    f[ic - ic/2 + i] = fi;
  }
  for (int i = -di; i <= di; i++)
  {
    const auto fi = max(ampl*(di - abs(i)),real_t{0});
    f[ic + ic/3 + i] = fi;
  }

}

template<typename vector_t>
void dump2file(const vector_t &f, const std::string fileName)
{
  const int n = f.size();
  std::ofstream fout(fileName);
  for (int i = 1; i < n-1; i++)
    fout << f[i] << std::endl;
}

int main(int argc, char * argv[])
{
  const int ncell = argc > 1 ? atoi(argv[1]) : 100;
  fprintf(std::cerr, "ncell= %\n", ncell);

  const int niter = argc > 2 ? atoi(argv[2]) : 10;
  fprintf(std::cerr, "niter= %\n", niter);

  const int ngs_iter = argc > 3 ? atoi(argv[3]) : 10;
  fprintf(std::cerr, "ngs_iter= %\n", ngs_iter);

  

  using real_t = double;

  using params_t = Params<real_t>;
  using vector_t = std::vector<real_t>;

  auto params = params_t{};
  params.dx   = 1;
  params.diff = 1;
  params.cfl  = 0.40;  /* stable for cfl <= 0.5 */
  params.cfl  = 8.0;

  vector_t f(ncell+2);

  set_ic(params,f);
  dump2file(f,"ic.txt");

  for (int iter = 0; iter < niter; iter++)
  {
    fprintf(std::cerr, "iter= %\n", iter);
#if 0
    compute_update<FreeBC>(params,f);
//    compute_update<PeriodicBC>(params,f);
#elif 1
    compute_update_cheby<FreeBC>(params,f,ngs_iter);
#endif
    dump2file(f,"iter"+std::to_string(iter));
  }




  return 0;  

}




