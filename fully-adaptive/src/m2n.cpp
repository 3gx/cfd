#include "LegendrePoly.h"

int main(int argc, char *argv[])
{
  using namespace std;

  constexpr int M = 2;
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
