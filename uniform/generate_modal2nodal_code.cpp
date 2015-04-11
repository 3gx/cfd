#include <iostream>
#include <string>
#include <cassert>
#include <array>
#include <cmath>
#include <algorithm>

template< bool B, class T = void >
using eIf = typename std::enable_if<B,T>::type;

template<typename... Ts>
struct IsEmptyBase
{
  bool value = sizeof...(Ts);
};

template<typename... Ts>
using IsEmpty = typename IsEmptyBase<Ts...>::value;

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
      return ((2*N-1)*(x*eval<N-1>(x)) - (N-1)*eval<N-2>(x))/N;
    }

  template<int P, typename T>
    static real_t evalR(T x)  
    { return eval<P>(static_cast<real_t>(x)); }

  template<int Pfirst, int... P, typename Tfirst, typename... T>
    static real_t eval(Tfirst x, T... args)
    {
      return evalR<Pfirst,Tfirst>(x)*eval<P...>(args...);
    }
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
    incr(F&& f)
    {
    }


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
    static void exec(F&& f)
    {
      warmup<DIM>(f);
    }
};

#if 0
template<int M, typename F>
struct static_loop3
{
  /*
   for (int a = 0; a <= M; a++)
    for (int b = 0; b <= M-a; b++)
      for (int c = 0; c <= M-a-b; c++)
      {
        ..
      }
   */
  template<int a, int b, int c,int COUNT>
    static eIf<(a>M)> eval(F &f)
    {
    }
  template<int a, int b, int c, int COUNT>
    static eIf<(a<=M && b>M-a)> eval(F &f)
    {
      eval<a+1,0,0,COUNT>(f);
    }
  template<int a, int b, int c, int COUNT>
    static eIf<(a<=M && b<= M-a && c>M-b-a)> eval(F &f)
    {
      eval<a,b+1,0,COUNT>(f);
    }

  template<int a, int b, int c, int COUNT = 0>
    static eIf<(a<=M && b<= M-a  && c<=M-b-a)> eval(F& f)
    {
      f. template eval<COUNT, c,b,a>();
      eval<a,b,c+1,COUNT+1>(f);
    }

  template<typename... T>
  static F loop(T... args)
  {
    F f(args...);
    eval<0,0,0>(f);
    return f;
  }
};
#endif

template<int M, typename real_t>
struct GenerateMatrix
{
  static constexpr int factorial(int n)
  {
    return n > 0 ? n*factorial(n-1) : 1;
  }
  static constexpr int product (int m, int i, int d)
  {
    return i < d ? (m+i+1)*product(m, i+1, d) : 1;
  }
  
  static constexpr int _size = product(M,0,3)/factorial(3);
  std::array<real_t,_size> matrix[_size];

  constexpr int size() const {return _size;}
  const real_t* getMatrix() const {return &matrix[0][0];}

  struct Expansion
  {
    std::array<real_t,_size> result;
    real_t x,y,z;
    Expansion(real_t _x, real_t _y, real_t _z) : x(_x), y(_y), z(_z) {};

    constexpr int size() const {return _size;}
    real_t operator[](const int i) const {return result[i];}

    template<int idx, int c,int b,int a>
      void eval()
      {
//        fprintf(stderr, "idx= %d M= %d : a= %d b= %d c= %d\n", idx, M, a,b,c);
        static_assert(idx < _size, "Buffer overflow");
        result[idx] = LegendrePoly<real_t>::template eval<a,b,c>(x,y,z);
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


    int count = 0;
    for (int a = 0; a <= M; a++)
      for (int b = 0; b <= M-a; b++)
        for (int c = 0; c <= M-a-b; c++)
        {
          assert(count < _size);
          const real_t x = nodes[c];
          const real_t y = nodes[b];
          const real_t z = nodes[a];
          Expansion f(x,y,z);
          static_loop<M,3>::exec(f);
          matrix[count] = f.result;
          count++;
        }
  }

  void printMatrix() const
  {
    for (int j = 0; j < _size; j++)
    {
      for (int i = 0; i < _size; i++)
      {
        printf("%5.2f ", matrix[j][i]);
      }
      printf("\n");
    }
  }
  void printMatrix(const real_t *matrix) const
  {
    for (int j = 0; j < _size; j++)
    {
      for (int i = 0; i < _size; i++)
      {
        printf("%5.2f ", matrix[j*_size + i]);
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
  dgetri_(&NN,Ainv,&NN,IPIV,WORK,&LLWORK,&INFO);
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


int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;
  constexpr int M = 4;
  using real_t = double;

#if 1  
  static_loop<M,3>::exec(GenerateMatrix<M,real_t>::Expansion(0.3,0.4,0.5));
  cout << "---------------\n";
  static_loop<M,4>::exec(Foo<M>{});

//  cout << product(4,0,3)/factorial(3) << endl;

#endif
  GenerateMatrix<M,real_t> g;


//  g.printMatrix(g.getMatrix());

  const int non_zero =
    std::count_if(g.getMatrix(), g.getMatrix()+g.size()*g.size(), [](const real_t val) {return std::abs(val) > 1.0e-10;});
  fprintf(stderr, " matrix size= [ %d x %d ]\n", g.size(), g.size());
  fprintf(stderr, " number of non-zero elements= %d [ %g %c ]\n", non_zero, non_zero*100.0/(g.size()*g.size()), '%' );

  real_t Ainv[g._size*g._size];
  inverse<g._size>(g.getMatrix(), Ainv);
  
//  g.printMatrix(Ainv);
  const int non_zero_inv =
    std::count_if(Ainv, Ainv+g.size()*g.size(), [](const real_t val) {return std::abs(val) > 1.0e-10;});
  fprintf(stderr, " number of non-zero-inv elements= %d [ %g %c ]\n", non_zero_inv, non_zero_inv*100.0/(g.size()*g.size()), '%' );

  fprintf(stderr, " -- tota number of non-zeros= %d [ %g %c ] \n",
      non_zero + non_zero_inv, (non_zero + non_zero_inv)*50.0/(g.size()*g.size()), '%');


  if (verify(g.getMatrix(), Ainv, g.size()))
    fprintf(stderr , " -- Matrix verified -- \n");
  cout << LegendrePoly<real_t>::template eval<1,2,3>(1.0,2.0,3.0) << endl;
  cout << LegendrePoly<real_t>::template eval<2>(0.0) << endl;

  const auto m2n = generateMatmulCode(g.getMatrix(), g.size());

  cout << m2n << endl;
  
  return 0;
}
