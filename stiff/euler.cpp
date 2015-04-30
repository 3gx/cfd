#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include "myprintf.h"

template<typename real_type>
auto solve(real_type h, real_type tmax, real_type y0)
{
  using pair_type = std::tuple<real_type,real_type>;
  using vector_type = std::vector<pair_type>;
  using std::get;

  const auto nstep = static_cast<int>(ceil(tmax/h));

  vector_type solution;

  auto time = 0.0;
  solution.push_back(pair_type{time,y0});

  auto rhs = [](const auto x) { return -15.0*x; };

  const auto scale= 1;
  const auto dt = h/scale;

  for (int i = 0; i < nstep*scale; i++)
  {
    const auto& cur = solution.back();
    const auto dy = rhs(get<1>(cur))*dt;
    const auto y1 = get<1>(cur) + dy;
    solution.emplace_back(get<0>(cur)+dt,y1);
  }


  return solution;

}

int main(int argc, char * argv[])
{
  using std::cout;
  using std::cerr;
  using std::endl;
  using std::get;
  const auto ncell = argc > 1 ? atoi(argv[1]) : 10;

  fprintf(cerr, "ncell= %\n", ncell);

  const auto h = 1.0/16/16;

  const auto tmax = 1.0;

  const auto y0 = 1.0;
  const auto solution = solve(h,tmax, y0);

  for (const auto &x : solution)
  {
    printf(" %15.15g %15.15g \n", get<0>(x), get<1>(x));
  }


  return 0;

}
