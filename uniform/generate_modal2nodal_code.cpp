#include <iostream>
#include <string>
#include <cassert>
#include <array>
#include <cmath>

template<typename real_t>
struct LegendrePoly
{
  template<int N>
    static typename std::enable_if<(N==1),real_t>::type eval(const real_t x) { return x; }

  template<int N>
    static typename std::enable_if<(N==0),real_t>::type eval(const real_t x) { return 1; }

  template<int N>
    static typename std::enable_if<(N>1),real_t>::type eval(const real_t x)
    {
      return ((2*N-1)*(x*eval<N-1>(x)) - (N-1)*eval<N-2>(x))/N;
    }


  template<int P, int Q, int R>
    static real_t Poly3D(const real_t x, const real_t y, const real_t z)
    {
      return eval<P>(x)*eval<Q>(y)*eval<R>(z);
    }
};

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
    static typename std::enable_if<(a>M),void>::type eval(F &f)
    {
    }
  template<int a, int b, int c, int COUNT>
    static typename std::enable_if<(a<=M && b>M-a),void>::type eval(F &f)
    {
      eval<a+1,0,0,COUNT>(f);
    }
  template<int a, int b, int c, int COUNT>
    static typename std::enable_if<(a<=M && b<= M-a && c>M-b-a),void>::type eval(F &f)
    {
      eval<a,b+1,0,COUNT>(f);
    }

  template<int a, int b, int c, int COUNT = 0>
    static typename std::enable_if<(a<=M && b<= M-a  && c<=M-b-a),void>::type eval(F& f)
    {
      f. template eval<a,b,c,COUNT>();
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

template<int M, typename real_t>
struct GenerateMatrix
{
  static constexpr int _size = (M+1)*(M+2)*(M+3)/6;
  std::array<real_t,_size> matrix[_size];

  constexpr int size() const {return _size;}
  const real_t* getMatrix() const {return &matrix[0][0];}

  struct Expansion3D
  {
    std::array<real_t,_size> result;
    real_t x,y,z;
    Expansion3D(real_t _x, real_t _y, real_t _z) : x(_x), y(_y), z(_z) {};

    constexpr int size() const {return result.size();}
    real_t operator[](const int i) const {return result[i];}

    template<int a,int b,int c, int idx>
      void eval()
      {
        static_assert(idx < _size, "Buffer overflow");
        result[idx] = LegendrePoly<real_t>::template Poly3D<a,b,c>(x,y,z);
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
          const auto f = static_loop3<M,Expansion3D>::loop(x,y,z);
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

#if 0
  Main(int argc, char *argv[])
  {
    using real_t = double;
    using P = LegendrePoly<real_t>;
    using std::cout;
    using std::endl;

    assert(argc > 1);
    const real_t x = atof(argv[1]);

    constexpr int N = 7;

    cout << "LP(" << N << ", x=" << x << ")= " <<  P::eval<N>(x) << endl;
    cout << "Root " << P::root<N>() << endl;

    cout << endl;

    const real_t xval=0.5;
    const real_t yval=0.5;
    const real_t zval=0.5;
    constexpr int M = 3;
    using loop1 = static_loop3<M,Expansion3D<M,real_t>>;
    int count = 0;
    const auto f = loop1::loop(xval,yval,zval);

    const int size = f.size();
    for (int i = 0; i < size; i++)
      fprintf(stderr, " P[%2d]= %g \n", i,f[i]);
  }
#endif


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


int main(int argc, char *argv[])
{
  constexpr int M = 3;
  using real_t = double;
  GenerateMatrix<M,real_t> g;

  g.printMatrix(g.getMatrix());

  const int non_zero =
    std::count_if(g.getMatrix(), g.getMatrix()+g.size()*g.size(), [](const real_t val) {return std::abs(val) > 1.0e-10;});
  fprintf(stderr, " matrix size= [ %d x %d ]\n", g.size(), g.size());
  fprintf(stderr, " number of non-zero elements= %d [ %g %c ]\n", non_zero, non_zero*100.0/(g.size()*g.size()), '%' );

  real_t Ainv[g.size()*g.size()];
  inverse<g.size()>(g.getMatrix(), Ainv);
  
  g.printMatrix(Ainv);
  const int non_zero_inv =
    std::count_if(Ainv, Ainv+g.size()*g.size(), [](const real_t val) {return std::abs(val) > 1.0e-10;});
  fprintf(stderr, " number of non-zero-inv elements= %d [ %g %c ]\n", non_zero_inv, non_zero_inv*100.0/(g.size()*g.size()), '%' );

  return 0;
}
