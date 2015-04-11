#include <iostream>
#include <string>
#include <cassert>
#include <array>
#include <cmath>
#include <algorithm>
#include <utility>

/* 
 * Code to generate code to compute modal-to-nodal and nodal-to-modal transform
 * Requires: C++14, lapack to compute inverse
 */

template< bool B, class T = void >
using enableIf = typename std::enable_if<B,T>::type;

/***************************************
 * Unpack parameter list from an array *
 ***************************************/
// C++11
template<size_t N, size_t... S>
struct Unpack : Unpack<N-1,N-1,S...>{};
template<size_t... S>
struct Unpack<0,S...>
{
  template<typename F, typename X>
  static auto eval(F&& f, X&& x) -> decltype(f(x[S]...))
  {
    return f(x[S]...);
  }
};
// C++14
template<typename F, typename X, size_t... I>
static auto unpack(F&& f, X&& x, std::index_sequence<I...>)
{
  return f(x[I]...);
}

/**************************************/

/* Compute Legendre Polynomials Pn(x) */
template<typename real_t>
struct LegendrePoly
{
  template<int N>
    static auto eval(const real_t x) -> enableIf<(N==1),real_t> 
    { return x; }

  template<int N>
    static auto eval(const real_t x) -> enableIf<(N==0),real_t>
    { return 1; }

  template<int N>
    static auto eval(const real_t x) -> enableIf<(N>1),real_t>
    {
      return ((2*N-1)*(x*eval<N-1>(x)) - (N-1)*eval<N-2>(x))*(1.0/N);
    }


#if 0  /* this won't work with type deduction. need explicit specilaization */
  template<int Pfirst, int... P, typename Tfirst, typename... T>
    static auto eval(Tfirst x, T... args) -> enableIf<(sizeof...(P)>0),real_t>
    {
      return eval<Pfirst>(x)*eval<P...>(args...);
    }
  template<int Pfirst, typename Tfirst>
    static auto eval(Tfirst x) -> real_t
    {
      return eval<Pfirst>(x);
    }
#else
  template<int P1, int P2, int P3, int P4>
    static real_t eval(real_t x, real_t y, real_t z, real_t w)
    {
      return eval<P1>(x)*eval<P2>(y)*eval<P3>(z)*eval<P4>(w);
    }
  template<int P1, int P2, int P3>
    static real_t eval(real_t x, real_t y, real_t z)
    {
      return eval<P1>(x)*eval<P2>(y)*eval<P3>(z);
    }
  template<int P1, int P2>
    static real_t eval(real_t x, real_t y)
    {
      return eval<P1>(x)*eval<P2>(y);
    }
#endif
};

template<size_t M, size_t DIM>
struct static_loop
{
  /* template meta-program for  this type of loop
   size_t count = 0;
   for (size_t a = 0; a <= M; a++)
    for (size_t b = 0; b <= M-a; b++)
      for (size_t c = 0; c <= M-a-b; c++)
        ...
      {
        f(count,a,b,c,...);
        count++;
      }
   */

  /******************/
  /* helper methods */
  /******************/
  template<size_t... Vs> 
    static constexpr auto sum() -> enableIf<sizeof...(Vs)==0,size_t> { return 0; }
  template<size_t V, size_t... Vs> 
    static constexpr auto sum() -> size_t { return V + sum<Vs...>(); }

  /**************/
  /* basic loop */
  /**************/
  template<size_t COUNT, size_t B, size_t... As, typename F>
    static auto eval(F&& f) -> enableIf<(B<=M-sum<As...>())> 
    {
      f.template eval<COUNT, B, As...>();  /* call function */
      eval<COUNT+1,B+1,As...>(f);
    }
  template<size_t COUNT, size_t B, size_t... As, typename F>
    static auto eval(F&& f) -> enableIf<(B>M-sum<As...>())> 
    {
      incr<COUNT,As...>(f, std::index_sequence<0>());
    }

  /*************/
  /* increment */
  /*************/
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, std::index_sequence<I...>) -> enableIf<(B<M-sum<As...>())>
    {
      eval<COUNT,I...,B+1,As...>(f);
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, std::index_sequence<I...>) -> enableIf<(B>=M-sum<As...>()) && (sizeof...(As) > 0)> 
    {
      incr<COUNT,As...>(f, std::index_sequence<0,I...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, std::index_sequence<I...>) -> enableIf<(B>=M-sum<As...>()) && (sizeof...(As) == 0)> 
    {}

  /***********/
  /* warm-up */ 
  /***********/
  template<size_t D, size_t... ZEROs, typename F>
    static auto warmup(F&& f,std::index_sequence<ZEROs...>) -> enableIf<(D>0)> 
    {
      warmup<D-1>(f,std::index_sequence<0,ZEROs...>());
    } 
  template<size_t D, size_t... ZEROs, typename F>
    static auto warmup(F&& f, std::index_sequence<ZEROs...>) -> enableIf<D==0> 
    {
      eval<0,ZEROs...>(f);
    } 

  /***************
   * entry point *
   ***************/
  template<typename F>
    static F& exec(F&& f)
    {
      warmup<DIM-1>(f,std::index_sequence<0>());
      return f;
    }
};

template<int M, int DIM, typename real_t>
struct GenerateMatrix
{
  static constexpr int factorial(int n)
  {
    return n > 0 ? n*factorial(n-1) : 1;
  }
  static constexpr int product (int m, int d, int i = 0)
  {
    return i < d ? (m+i+1)*product(m, d, i+1) : 1;
  }
  
  static constexpr int size = product(M,DIM,0)/factorial(DIM);
  std::array<real_t,size> matrix[size];

  const real_t* getMatrix() const {return &matrix[0][0];}

  /* Helper to compute matrix raw for a given node */
  struct Expansion
  {
    std::array<real_t,DIM > node;
    std::array<real_t,size> result;

    Expansion(const std::array<real_t,DIM>& v) : node(v) { std::reverse(node.begin(), node.end()); }

    real_t operator[](const int i) const {return result[i];}

    template<int count, int... Vs>
      void eval()
      {
        static_assert(count < size, "Buffer overflow");
//        result[count] = Unpack<DIM>::eval(LegendrePoly<real_t>::template eval<Vs...>,node);
        result[count] = unpack(LegendrePoly<real_t>::template eval<Vs...>,node,std::make_index_sequence<DIM>());
      }
  };

  /* Help that computes a map between a linear index and (abc..) indices */
  struct Indices
  {
    using index_t = std::array<int,DIM>;
    std::array<index_t,size> indices;

    const index_t& operator[](const int i) const {return indices[i];}

    template<int count, int V, int... Vs>
      void fill()
      {
        indices[count][sizeof...(Vs)] = V;
        fill<count, Vs...>();
      }
    template<int count, int... Vs>
      auto fill() -> enableIf<sizeof...(Vs)==0>
      {}

    template<int count, int... Vs>
      void eval()
      {
        static_assert(count < size, "Buffer overflow");
        fill<count, Vs...>();
      }
  };

  /* Compute matrix transformation from modal to nodal expansion */
  GenerateMatrix()
  {
    real_t nodes[] = 
    {
      /* n11 */ 0.0,
      /* n21 */ -0.57735026918962576451,
      /* n22 */ +0.57735026918962576451,
      /* n31 */ -0.77459666924148337704,
      /* n33 */ +0.77459666924148337704,
      /* n41 */ -0.86113631159405257522,
      /* n44 */ +0.86113631159405257522,
      /* n42 */ -0.33998104358485626480,
      /* n43 */ +0.33998104358485626480,
      /* n51 */ -0.90617984593866399280,
      /* n51 */ +0.90617984593866399280
      /* n52 */ -0.53846931010568309104,
      /* n52 */ +0.53846931010568309104,
    };
    static_assert(M>=0 && M<13, "M-value is out of range");

    const auto indices = static_loop<M,DIM>::exec(Indices{});
    std::array<real_t,DIM> node;

    for (int i = 0; i < size; i++)
    {
      const auto& idx = indices[i];
      for (int k = 0; k < DIM; k++)
        node[k] = nodes[idx[DIM-k-1]];
      const auto f = static_loop<M,DIM>::exec(Expansion{node});
      matrix[i] = f.result;
    }
  }

  void printMatrix(const real_t *matrix) const
  {
    for (int j = 0; j < size; j++)
    {
      for (int i = 0; i < size; i++)
      {
        printf("%5.2f ", matrix[j*size + i]);
      }
      printf("\n");
    }
  }
  
  /* prints matrix */
  void printMatrix() const
  {
    printMatrix(getMatrix());
  }
};

template<typename real_t>
std::string generateMatmulCode(const real_t *matrix, const int size)
{
  std::string code;

  return code;
};

extern "C" 
{
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

template<int N>
void inverse(const double* A, double *Ainv)
{
  constexpr int LWORK = N*N;

  int IPIV[N+1];
  double WORK[LWORK];
  int INFO;

  for (int i = 0; i < N*N; i++)
    Ainv[i] = A[i];

  int NN = N;
  int LLWORK = LWORK;

  dgetrf_(&NN,&NN,Ainv,&NN,IPIV,&INFO);
  assert(INFO == 0);
  dgetri_(&NN,Ainv,&NN,IPIV,WORK,&LLWORK,&INFO);
  assert(INFO == 0);
}

bool verify(const double *A, const double *B, const int N)
{
  bool success = true;
  constexpr double eps = 1.0e-13;
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      double res = 0;
      for (int k = 0; k < N; k++)
        res += A[i*N+k]*B[k*N+j];
      if (
          ((i==j) && std::abs(res - 1.0) > eps) ||
          ((i!=j) && std::abs(res)       > eps)
         )
      {
        printf(" (%2d,%2d) = %5.2f \n", i,j, res);
        success = false;
      }
    }
  }
  return success;
}


template<int M>
struct Printer
{
  template<int A, int... As>
    auto evalR() -> enableIf<(sizeof...(As)>0)>
    {
      fprintf(stderr, "%c= %d ", static_cast<char>('a'+sizeof...(As)), A);
      evalR<As...>();
    }
  
  template<int A, int... As>
    auto evalR() -> enableIf<(sizeof...(As)==0)>
    {
      fprintf(stderr, "%c= %d \n", 'a', A);
    }

    template<int count, int... As>
      void eval()
      {
        fprintf(stderr, "idx= %d M= %d : ", count, M);
        evalR<As...>();
      }
};

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  constexpr int M = 4;
  constexpr int DIM = 4;
  using real_t = double;

  cout << LegendrePoly<real_t>::template eval<1,2,3>(1.0,2.0,3.0) << endl;
  cout << Unpack<3>::template eval(LegendrePoly<real_t>::template eval<1,2,3>,std::array<real_t,3>{1,2,3}) << endl;

//  return 0;

  static_loop<M,DIM>::exec(Printer<M>{});

  GenerateMatrix<M,DIM,real_t> g;


//  g.printMatrix(g.getMatrix());

  const int non_zero =
    std::count_if(g.getMatrix(), g.getMatrix()+g.size*g.size, [](const real_t val) {return std::abs(val) > 1.0e-10;});
  fprintf(stderr, " matrix size= [ %d x %d ]\n", g.size, g.size);
  fprintf(stderr, " number of non-zero elements= %d [ %g %c ]\n", non_zero, non_zero*100.0/(g.size*g.size), '%' );

  real_t Ainv[g.size*g.size];
  inverse<g.size>(g.getMatrix(), Ainv);
  
//  g.printMatrix(Ainv);
  const int non_zero_inv =
    std::count_if(Ainv, Ainv+g.size*g.size, [](const real_t val) {return std::abs(val) > 1.0e-10;});
  fprintf(stderr, " number of non-zero-inv elements= %d [ %g %c ]\n", non_zero_inv, non_zero_inv*100.0/(g.size*g.size), '%' );

  fprintf(stderr, " -- total number of non-zeros= %d [ %g %c ] \n",
      non_zero + non_zero_inv, (non_zero + non_zero_inv)*50.0/(g.size*g.size), '%');


  if (verify(g.getMatrix(), Ainv, g.size))
    fprintf(stderr , " -- Matrix verified -- \n");
//  cout << LegendrePoly<real_t>::template eval<1,2,3>(1.0,2.0,3.0) << endl;
//  cout << LegendrePoly<real_t>::template eval<2>(0.0) << endl;
//
  fprintf(stderr, " M= %d  DIM= %d \n", M, DIM);

  const auto m2n = generateMatmulCode(g.getMatrix(), g.size);

  cout << m2n << endl;
  
  return 0;
}
