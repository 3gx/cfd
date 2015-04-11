#include <iostream>
#include <string>
#include <cassert>
#include <array>
#include <cmath>
#include <algorithm>
#include <utility>

template< bool B, class T = void >
using eIf = typename std::enable_if<B,T>::type;


template<typename real_t>
struct LegendrePoly
{
  template<int N>
    static eIf<(N==1),real_t> eval(const real_t x) { return x; }

  template<int N>
    static eIf<(N==0),real_t> eval(const real_t x) { return 1; }

  template<int N>
    static eIf<(N>1),real_t> eval(const real_t x)
    {
      return ((2*N-1)*(x*eval<N-1>(x)) - (N-1)*eval<N-2>(x))*(1.0/N);
    }

  template<int P, typename T>
    static real_t evalR(T x)  
    { return eval<P>(static_cast<real_t>(x)); }
  template<int Pfirst, int... P, typename Tfirst, typename... T>
    static real_t eval(Tfirst x, T... args)
    {
      return evalR<Pfirst,Tfirst>(x)*eval<P...>(args...);
    }

#if 1
  template<int Pfirst, int... Ps>
    static real_t eval(const std::initializer_list<real_t>& arg)
    {
      return eval<Pfirst>(*(arg.end() - sizeof...(Ps)-1))*eval<Ps...>(arg);
    }

  template<int... Ps>
    static auto eval(const std::initializer_list<real_t>& arg) -> eIf<sizeof...(Ps)==0,real_t>
    { return 1.0;}
#endif

#if 1
  template<int Pfirst, int... Ps, typename Arg>
    static real_t eval(Arg&& arg)
    {
      return eval<Pfirst>(*(arg.end() - sizeof...(Ps)-1))*eval<Ps...>(std::forward<Arg>(arg));
    }

  template<int... Ps, typename Arg>
    static auto eval(Arg&& arg) -> eIf<sizeof...(Ps)==0,real_t>
    { return 1.0;}
#endif
};

template<int M, int DIM>
struct static_loop
{
  /* template meta-program for  this type of loop
   int count = 0;
   for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M-a; b++)
      for (int c = 0; c <= M-a-b; c++)
        ...
      {
        f(count,c,b,a);
        count++;
      }
   */

  /******************/
  /* helper methods */
  /******************/
  template<int... Vs> 
    static constexpr eIf<sizeof...(Vs)==0,int> sum()  { return 0; }
  template<int V, int... Vs> 
    static constexpr int sum() { return V + sum<Vs...>(); }

  /**************/
  /* basic loop */
  /**************/
  template<int COUNT, int B, int... As, typename F>
    static eIf<(B<=M-sum<As...>())> eval(F&& f)
    {
      f.template eval<COUNT, B, As...>();  /* call function */
      eval<COUNT+1,B+1,As...>(f);
    }
  template<int COUNT, int B, int... As, typename F>
    static eIf<(B>M-sum<As...>())> eval(F&& f)
    {
      incr<1, COUNT, As...>(f);
    }

  /*************/
  /* increment */
  /*************/
  template<int K, int COUNT, int B, int... As, typename F>
    static eIf<(B<M-sum<As...>())> incr(F&& f)
    {
      cont<K,COUNT,B+1,As...>(f);
    }
  template<int K, int COUNT, int B, int... As, typename F>
    static eIf<(B>=M-sum<As...>()) && (sizeof...(As) > 0)> 
    incr(F&& f)
    {
      incr<K+1,COUNT,As...>(f);
    }
  template<int K, int COUNT, int B, int... As, typename F>
    static eIf<(B>=M-sum<As...>()) && (sizeof...(As) == 0)> 
    incr(F&& f) {}


  /************/
  /* continue */
  /************/
  template<int K, int COUNT, int... As, typename F>
    static eIf<(K>0)> cont(F&& f)
    {
      cont<K-1,COUNT,0,As...>(f);
    }
  template<int K, int COUNT, int... As, typename F>
    static eIf<K==0> cont(F&& f)
    {
      eval<COUNT, As...>(f);
    }

  /***********/
  /* warm-up */ 
  /***********/
  template<int D, int... ZERO, typename F>
    static eIf<(D>0)> warmup(F&& f)
    {
      warmup<D-1,0,ZERO...>(f);
    } 
  template<int D, int... ZERO, typename F>
    static eIf<D==0> warmup(F&& f)
    {
      eval<D,ZERO...>(f);
    } 

  /***************
   * entry point *
   ***************/
  template<typename F>
    static F& exec(F&& f)
    {
      warmup<DIM>(f);
      return f;
    }
};

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

template<int M, int DIM, typename real_t>
struct GenerateMatrix
{
//  static constexpr int DIM = 3;
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
        result[count] = LegendrePoly<real_t>::template eval<Vs...>(node); 
//        Unpack<DIM>::template eval(LegendrePoly<real_t>::template eval<Vs...>,node);
      }
  };

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
      eIf<sizeof...(Vs)==0> fill() {}

    template<int count, int... Vs>
      void eval()
      {
        static_assert(count < size, "Buffer overflow");
        fill<count, Vs...>();
      }

  };

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
      /* n42 */ -0.33998104358485626480,
      /* n43 */ +0.33998104358485626480,
      /* n44 */ +0.86113631159405257522,
      /* n51 */ -0.90617984593866399280,
      /* n52 */ -0.53846931010568309104,
      /* n52 */ +0.53846931010568309104,
      /* n51 */ +0.90617984593866399280
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

  void printMatrix() const
  {
    for (int j = 0; j < size; j++)
    {
      for (int i = 0; i < size; i++)
      {
        printf("%5.2f ", matrix[j][i]);
      }
      printf("\n");
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
struct Foo
{
  template<int A, int... As>
    eIf<(sizeof...(As)>0)>
    evalR()
    {
      fprintf(stderr, "%c= %d ", static_cast<char>('a'+sizeof...(As)), A);
      evalR<As...>();
    }
  
  template<int A, int... As>
    eIf<(sizeof...(As)==0)>
    evalR()
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


template<typename V>
struct Car
{
  template<int W>
    static V foo(int x, int y, int z, int v, int w)
    {
      using std::cout;
      using std::endl;
      cout << x << endl;
      cout << y << endl;
      cout << z << endl;
      cout << v << endl;
      cout << w << endl;
      return V(100.3)+W+x+y+z+w+v+w;
    }
};

template<typename V, int W>
struct Mar
{
  static void foo()
  {
    std::array<int,5> var = {1,2,3,4,5};
    std::cout << Unpack<5>::template eval(Car<V>::template foo<W>,var) << std::endl;
  }
};

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  constexpr int M = 3;
  constexpr int DIM = 4;
  using real_t = double;

  std::array<int,5> var = {1,2,3,4,5};
  auto f = [&](int x, int y, int z, int v, int w) -> int
  {
    cout << x << endl;
    cout << y << endl;
    cout << z << endl;
    cout << v << endl;
    cout << w << endl;
    return x+y+z+w+v+w;
  };
  Mar<float,1000>::foo();
  cout << LegendrePoly<real_t>::template eval<1,2,3>(1.0,2.0,3.0) << endl;
  cout << LegendrePoly<real_t>::template eval<1,2,3>({1.0,2.0,3.0}) << endl;
  


//  return 0;

#if 1  
  static_loop<M,4>::exec(Foo<M>{});

//  cout << product(4,0,3)/factorial(3) << endl;

#endif
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
