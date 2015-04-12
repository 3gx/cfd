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

/* helpers */
template< bool B, class T = void >
using enableIf = typename std::enable_if<B,T>::type;

template<size_t... I>
using indexSeq = std::index_sequence<I...>;

template<size_t N>
using makeIndexSeq = std::make_index_sequence<N>;


/***************************************
 * Unpack parameter list from an array *
 ***************************************/
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

  static constexpr real_t getRoot(const int n)
  {
#if 0 /* non-return constexpr are not supported yet in GCC <= 4.9.2 */
      constexpr real_t roots[] = 
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

      return roots[n];
#else
      return 
        n ==  0 ? /* n11 */ 0.0 :
        n ==  1 ? /* n21 */ -0.57735026918962576451 :
        n ==  2 ? /* n22 */ +0.57735026918962576451 :
        n ==  3 ? /* n31 */ -0.77459666924148337704 :
        n ==  4 ? /* n33 */ +0.77459666924148337704 :
        n ==  5 ? /* n41 */ -0.86113631159405257522 :
        n ==  6 ? /* n44 */ +0.86113631159405257522 :
        n ==  7 ? /* n42 */ -0.33998104358485626480 :
        n ==  8 ? /* n43 */ +0.33998104358485626480 :
        n ==  9 ? /* n51 */ -0.90617984593866399280 :
        n == 10 ? /* n51 */ +0.90617984593866399280 :
        n == 11 ? /* n52 */ -0.53846931010568309104 :
        n == 12 ? /* n52 */ +0.53846931010568309104 :
        -1.0e38;
#endif
  }
};

template<size_t M, size_t DIM, typename Indices = makeIndexSeq<M>>
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
      f.template eval<COUNT, B, As...>(Indices());  /* call function */
      eval<COUNT+1,B+1,As...>(f);
    }
  template<size_t COUNT, size_t B, size_t... As, typename F>
    static auto eval(F&& f) -> enableIf<(B>M-sum<As...>())> 
    {
      incr<COUNT,As...>(f, indexSeq<0>());
    }

  /*************/
  /* increment */
  /*************/
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B<M-sum<As...>())>
    {
      eval<COUNT,I...,B+1,As...>(f);
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B>=M-sum<As...>()) && (sizeof...(As) > 0)> 
    {
      incr<COUNT,As...>(f, indexSeq<0,I...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B>=M-sum<As...>()) && (sizeof...(As) == 0)> 
    {}

  /***********/
  /* warm-up */ 
  /***********/
  template<size_t D, size_t... ZEROs, typename F>
    static auto warmup(F&& f,indexSeq<ZEROs...>) -> enableIf<(D>0)> 
    {
      warmup<D-1>(f,indexSeq<0,ZEROs...>());
    } 
  template<size_t D, size_t... ZEROs, typename F>
    static auto warmup(F&& f, indexSeq<ZEROs...>) -> enableIf<D==0> 
    {
      eval<0,ZEROs...>(f);
    } 

  /***************
   * entry point *
   ***************/
  template<typename F>
    static F& exec(F&& f)
    {
      warmup<DIM-1>(f,indexSeq<0>());
      return f;
    }
};

template<int M, int DIM, typename real_t, typename Poly = LegendrePoly<real_t>>
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

#if 0
    Expansion(const std::array<real_t,DIM>& v) : node(v) { std::reverse(node.begin(), node.end()); }
#else
    template<typename... Ts>
      Expansion(Ts... ts) : node{ts...} { std::reverse(node.begin(), node.end()); }
#endif

    real_t operator[](const int i) const {return result[i];}

    template<int count, int... Vs, size_t... I>
      void eval(indexSeq<I...>)
      {
        static_assert(count < size, "Buffer overflow");
        result[count] = Poly::template eval<Vs...>(node[I]...);
      }
  };

  /* Help that computes a map between a linear index and (abc..) indices */
  struct Indices
  {
    using index_t = std::array<int,DIM>;
    std::array<index_t,size> indices;

    const index_t& operator[](const int i) const {return indices[i];}

    template<int count, int... Vs, size_t... I>
      void eval(indexSeq<I...>)
      {
        static_assert(count < size, "Buffer overflow");
        auto tmp = {(indices[count][sizeof...(I)-1-I] = Vs)...};
      }
  };

  /* Compute matrix transformation from modal to nodal expansion */
  GenerateMatrix()
  {
    fillMatrixWithRoots(makeIndexSeq<DIM>());
  }

  template<size_t... I>
    void fillMatrixWithRoots(indexSeq<I...>)
    {
      const auto indices = static_loop<M,DIM>::exec(Indices{});
      for (int i = 0; i < size; i++)
      {
        const auto& idx = indices[i];
        const auto f = static_loop<M,DIM>::exec(Expansion(Poly::getRoot(idx[DIM-1-I])...));
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
        printf(" %sI[%2d,%2d] = %g \n", (i==j ? "1-" : "  "), i,j, 
            (i==j ? 1.0-res : res));
        success = false;
      }
    }
  }
  return success;
}


struct Printer
{
  template<int count, int... As, size_t... I>
    void eval(indexSeq<I...>)
    {
      const auto map = {std::tuple<char,int>{static_cast<char>('a'+sizeof...(I)-1-I),static_cast<int>(As)}...};

      fprintf(stderr, "idx= %d M= %d : ", count, (int)sizeof...(As));
      for (const auto& m : map)
      {
        fprintf(stderr, "%c= %d ", std::get<0>(m), std::get<1>(m));
      }
      fprintf(stderr, " \n");
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

  static_loop<M,DIM>::exec(Printer{});

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
