#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
#include <array>
#include <cmath>
#include <algorithm>
#include <vector>
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

/* genrate Zero index sequence */
template<size_t N, size_t... S>
struct makeZeroIndexSeqImpl : makeZeroIndexSeqImpl<N-1,0,S...> {};
template<size_t... S>
struct makeZeroIndexSeqImpl<0,S...>
{
  using type = indexSeq<S...>;
};
template<size_t N>
using makeZeroIndexSeq = typename makeZeroIndexSeqImpl<N>::type;

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


#if 0  /* self-reminder:
          this won't work with type deduction. need explicit specilaization */
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
  template<int P1, int P2, int P3, int P4, int P5>
    static real_t eval(real_t x, real_t y, real_t z, real_t w, real_t v)
    {
      return eval<P1>(x)*eval<P2>(y)*eval<P3>(z)*eval<P4>(w)*eval<P5>(v);
    }
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
#if 0 /* self-reminder:
         non-return constexpr are not supported yet in GCC <= 4.9.2 */
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

template<typename real_t, size_t N>
struct Vector
{
  using value_type = real_t;
  std::array<real_t,N> data;

  real_t& operator[](const size_t i) {return data[i];}
  const real_t& operator[](const size_t i) const {return data[i];}

  constexpr size_t size() const { return N; }
  
  auto begin ()       {return data.begin();}
  auto end   ()       {return data.end  ();}
  auto begin () const {return data.begin();}
  auto end   () const {return data.end  ();}
  auto cbegin() const {return data.cbegin();}
  auto cend  () const {return data.cend  ();}
  
  auto&&       front()       {return data.front();}
  const auto&& front() const {return data.front();}
  auto&&       back ()       {return data.back ();}
  const auto&& back () const {return data.back ();}

};

template<typename real_t, size_t M, size_t N=M, typename vector_t = Vector<real_t,N>>
struct Matrix
{
  using value_type = real_t;
  using vector_type = vector_t;

  std::array<vector_t,M> data;
  
  constexpr size_t        size() const {return M;}
  constexpr size_t vector_size() const {return N;}
  
  vector_t& operator[](const size_t i) {return data[i];}
  const vector_t& operator[](const size_t i) const {return data[i];}

  auto begin ()       {return data.front().begin();}
  auto end   ()       {return data.back ().end  ();}
  auto begin () const {return data.front().begin();}
  auto cend  () const {return data.back ().end  ();}
  auto cbegin() const {return data.front().begin();}
  auto end   () const {return data.back ().end  ();}
};



template<size_t M, size_t DIM, typename Indices = makeIndexSeq<DIM>>
struct UniqueStaticLoop
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

  static constexpr int factorial(int n)
  {
    return n > 0 ? n*factorial(n-1) : 1;
  }
  static constexpr int product (int m, int d, int i = 0)
  {
    return i < d ? (m+i+1)*product(m, d, i+1) : 1;
  }
  static constexpr int size = product(M,DIM,0)/factorial(DIM);

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
  template<size_t COUNT = 0, size_t B, size_t... As, typename F>
    static auto eval(F&& f, indexSeq<B,As...>) -> enableIf<(B<=M-sum<As...>())> 
    {
      f.template eval<COUNT, B, As...>(Indices());  /* call function */
      eval<COUNT+1>(f,indexSeq<B+1,As...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F>
    static auto eval(F&& f,indexSeq<B,As...>) -> enableIf<(B>M-sum<As...>())> 
    {
      incr<COUNT,As...>(f, indexSeq<0>());
    }

  /*************/
  /* increment */
  /*************/
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B<M-sum<As...>())>
    {
      eval<COUNT>(f,indexSeq<I...,B+1,As...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B>=M-sum<As...>()) && (sizeof...(As) > 0)> 
    {
      incr<COUNT,As...>(f, indexSeq<0,I...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B>=M-sum<As...>()) && (sizeof...(As) == 0)> 
    {}

  /***************
   * entry point *
   ***************/
  template<typename F>
    static F& exec(F&& f)
    {
      eval(f,makeZeroIndexSeq<DIM>());
      return f;
    }
};

template<size_t M, size_t DIM, typename Indices = makeIndexSeq<DIM>>
struct FullStaticLoop
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

  static constexpr int product (int m, int d)
  {
    return d > 0 ? (m+1)*product(m, d-1) : 1;
  }
  static constexpr int size = product(M,DIM);

  /**************/
  /* basic loop */
  /**************/
  template<size_t COUNT = 0, size_t B, size_t... As, typename F>
    static auto eval(F&& f, indexSeq<B,As...>) -> enableIf<(B<=M)> 
    {
      f.template eval<COUNT, B, As...>(Indices());  /* call function */
      eval<COUNT+1>(f,indexSeq<B+1,As...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F>
    static auto eval(F&& f,indexSeq<B,As...>) -> enableIf<(B>M)> 
    {
      incr<COUNT,As...>(f, indexSeq<0>());
    }

  /*************/
  /* increment */
  /*************/
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B<M)>
    {
      eval<COUNT>(f, indexSeq<I...,B+1,As...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B>=M) && (sizeof...(As) > 0)> 
    {
      incr<COUNT,As...>(f, indexSeq<0,I...>());
    }
  template<size_t COUNT, size_t B, size_t... As, typename F, size_t... I>
    static auto incr(F&& f, indexSeq<I...>) -> enableIf<(B>=M) && (sizeof...(As) == 0)> 
    {}

  /***************
   * entry point *
   ***************/
  template<typename F>
    static F& exec(F&& f)
    {
      eval(f,makeZeroIndexSeq<DIM>());
      return f;
    }
};

template<int M, int DIM, typename real_t, typename Poly = LegendrePoly<real_t>, typename StaticLoop = UniqueStaticLoop<M,DIM>, typename matrix_t = Matrix<real_t,StaticLoop::size>>
struct GenerateMatrix
{
  static constexpr int size = StaticLoop::size; 

  using static_loop_type = StaticLoop; 
  using matrix_type = matrix_t;
  using vector_t = typename matrix_t::vector_type;

  matrix_t matrix;

  const matrix_t& getMatrix() const {return matrix;}

  /* Helper to compute matrix raw for a given node */
  struct Expansion
  {
    std::array<real_t,DIM > node;
    vector_t result;

    template<typename... Ts>
      Expansion(Ts... ts) : node{ts...} { std::reverse(node.begin(), node.end()); }

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
      auto eval(indexSeq<I...>)
      {
        static_assert(count < size, "Buffer overflow");
        auto tmp = {(indices[count][sizeof...(I)-1-I] = Vs)...};
        return tmp;
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
      const auto indices = StaticLoop::exec(Indices{});
      for (int i = 0; i < size; i++)
      {
        const auto& idx = indices[i];
        const auto f = StaticLoop::exec(Expansion(Poly::getRoot(idx[DIM-1-I])...));
        matrix[i] = f.result;
      }
    }
};


template<typename matrix_t>
void printMatrix(const matrix_t& matrix) 
{
  const int size = matrix.size();
  for (int j = 0; j < size; j++)
  {
    for (int i = 0; i < size; i++)
    {
      printf("%5.2f ", matrix[j][i]);
    }
    printf("\n");
  }
}

template<typename matrix_t, typename real_t = typename matrix_t::value_type>
size_t countNonZeros(const matrix_t& matrix, const real_t eps = 1.0e-12)
{
  return std::count_if(
      matrix.begin(),
      matrix.end(),
      [eps](const real_t val) {return std::abs(val) > eps;});
}

template<typename matrix_t, typename real_t = typename matrix_t::value_type>
bool verifyMatrix(const matrix_t& A, const matrix_t& B, const real_t eps = 1.0e-12)
{
  const int N = A.size();
  bool success = true;
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      double res = 0;
      for (int k = 0; k < N; k++)
        res += A[i][k]*B[k][j];
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

extern "C" 
{
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

template<typename matrix_t, typename real_t = typename matrix_t::value_type>
auto invertMatrix(const matrix_t A) -> enableIf<std::is_same<real_t,double>::value, matrix_t>
{
  matrix_t Ainv;
  constexpr int N = A.size();
  constexpr int LWORK = N*N;

  int IPIV[N+1];
  double WORK[LWORK];
  int INFO;

  std::copy(A.begin(), A.end(), Ainv.begin());

  int NN = N;
  int LLWORK = LWORK;

  dgetrf_(&NN,&NN,Ainv.begin(),&NN,IPIV,&INFO);
  assert(INFO == 0);
  dgetri_(&NN,Ainv.begin(),&NN,IPIV,WORK,&LLWORK,&INFO);
  assert(INFO == 0);

  return Ainv;
}


template<typename matrix_t, typename real_t = typename matrix_t::value_type>
std::string generateMatvecCode(const matrix_t matrix, const real_t eps = 1.0e-12)
{
  using namespace std;
  ostringstream code;

  vector<string> bvar, xvar;

  const int size = matrix.size();
  for (int i = 0; i < size; i++)
  {
    const auto num = std::to_string(i);
    bvar.push_back("b["+num+"]");
    xvar.push_back("x["+num+"]");
  }

  for (int j = 0; j < size; j++)
  {
    code << xvar[j] << " = ";

    for (int i = 0; i < size; i++)
    {
      const auto value = matrix[j][i];
      char buf[256] = {0};
      if (std::abs(value) > eps)
      {
        sprintf(buf," + (%a)*%s",value,bvar[i].c_str());
        code << string(buf);
      }
    }
    code << ";" << endl;
  }

  return code.str();
};

struct Printer
{
  template<int count, int... As, size_t... I>
    void eval(indexSeq<I...>)
    {
      fprintf(stderr, "idx= %d M= %d : ", count, (int)sizeof...(As));

#ifdef __GNUC__  /* bug in GCC <= 4.9.2 tuple ? */
      const auto map = {(std::pair<char,int>{static_cast<char>('a'+sizeof...(I)-1-I),static_cast<int>(As)})...};
#else  /* __clang__ */
      const auto map = {(std::tuple<char,int>{static_cast<char>('a'+sizeof...(I)-1-I),static_cast<int>(As)})...};
#endif
      for (const auto& m : map)
      {
        fprintf(stderr, "%c= %d ", std::get<0>(m), std::get<1>(m));
      }

      fprintf(stderr, " \n");
    }
};
int main(int argc, char *argv[])
{
  using namespace std;

  constexpr int M = 4;
  constexpr int DIM = 4;
  using real_t = double;

#ifndef FULLLOOP
  using generate_matrix_t = GenerateMatrix<M,DIM,real_t>;
#else
  using generate_matrix_t = GenerateMatrix<M,DIM,real_t,LegendrePoly<real_t>,FullStaticLoop<M,DIM>>;
#endif
  generate_matrix_t g;
  generate_matrix_t::static_loop_type::exec(Printer{});

  if (argc > 1)
    printMatrix(g.getMatrix());
  
  /*** relative accuracy paramter ***/
  const real_t eps = 1.0e-12;      

  const int non_zero = countNonZeros(g.getMatrix(), eps);
  fprintf(stderr, " matrix size= [ %d x %d ]\n", g.size, g.size);
  fprintf(stderr, " number of non-zero elements= %d [ %g %c ]\n", non_zero, non_zero*100.0/(g.size*g.size), '%' );

  const auto Ainv = invertMatrix(g.getMatrix());
 
  if (argc > 2) 
    printMatrix(Ainv);
  const int non_zero_inv =countNonZeros(Ainv, eps);
  fprintf(stderr, " number of non-zero-inv elements= %d [ %g %c ]\n", non_zero_inv, non_zero_inv*100.0/(g.size*g.size), '%' );
  fprintf(stderr, " -- total number of non-zeros= %d [ %g %c ] \n",
      non_zero + non_zero_inv, (non_zero + non_zero_inv)*50.0/(g.size*g.size), '%');



  /*********************
   *** verify matrix ***
   *********************/

  if (verifyMatrix(g.getMatrix(), Ainv,eps))
    fprintf(stderr , " -- Matrix verified -- \n");
  else
  {
    fprintf(stderr, " -- Verification failure: exit -- \n");
    return -1;
  }
  fprintf(stderr, " M= %d  DIM= %d \n", M, DIM);


  /*********************************************
   ***** generate and write code to a file *****
   *********************************************/

  auto writeCode = [eps](auto matrix, auto fileName)
  {
    std::ofstream fout(fileName);
    fout << generateMatvecCode(matrix,eps) << std::endl;
  };
  writeCode(g.getMatrix(), "m2n.h");
  writeCode(Ainv,          "n2m.h");

  return 0;
}
