#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <array>
#include "myprintf.h"


template<typename T, size_t N>
constexpr auto array_size(T(&)[N]) { return N; }


template<typename real_type>
auto solve3(real_type h, real_type tmax, real_type y0)
{
  using pair_type = std::tuple<real_type,real_type>;
  using vector_type = std::vector<pair_type>;
  using std::get;

  const auto nstep = static_cast<int>(ceil(tmax/h));

  vector_type solution;
  
  real_type weights[] = {0.2777777777778,  0.4444444444444, 0.27777777777778};

  auto computeSolvec = [](auto C, auto h) 
  {
    const auto det = 0.0342935528120713*(60.0000000000000+C*h*(36.0000000000000+C*h*(9.00000000000000+C*h)));
    const auto idet = 1.0/det;
    return std::array<decltype(C),3>{
      idet*(0.00411522633744856*(500.000000000000+5.00000000000000*C*h*(48.7298334620742+8.87298334620742*C*h))),
      idet*(-0.0514403292181070*(-40.0000000000000+C*h*(-4.00000000000000+C*h))),
      idet*(0.00411522633744856*(500.000000000000-5.00000000000000*C*h*(28.7298334620742-1.12701665379258*C*h)))
    };
  };

  auto time = 0.0;
  solution.push_back(pair_type{time,y0});

  const auto scale = 1;
  const auto dt = h/scale;

  const auto C = 15.0;
  for (int i = 0; i < nstep*scale; i++)
  {
    const auto& u0 = solution.back();

    const auto solvec = computeSolvec(C,dt);
    auto dy = 0.0;
    for (int i = 0; i < (int)solvec.size(); i++)
    {
      const auto coeff = solvec[i] * get<1>(u0) * dt;
      dy += weights[i]*coeff;
    }
    dy *= -C;

    solution.emplace_back(get<0>(u0)+dt,get<1>(u0)+dy);
  }


  return solution;

}

template<typename real_type>
auto solve5(real_type h, real_type tmax, real_type y0)
{
  using pair_type = std::tuple<real_type,real_type>;
  using vector_type = std::vector<pair_type>;
  using std::get;

  const auto nstep = static_cast<int>(ceil(tmax/h));

  vector_type solution;
  
  real_type weights[] ={0.118463442528,0.239314335250,0.284444444444,0.239314335250,0.118463442528};

  auto computeSolvec = [](auto C, auto h) 
  {
    const auto det = 
      15120. + C*h*(8400. + C*h*(2100. + C*h*(300. + C*h*(25. + 1.*C*h))));
    const auto idet = 1.0/det;
    return std::array<decltype(C),5>{
      idet*(15120. + C*h*(7690.72 + C*h*(1722.59 + C*h*(210.471 + 13.0961*C*h)))),
      idet*(15120. + C*h*(4910.83 + C*h*(564.161 + C*h*(8.08594 - 3.73216*C*h)))),
      idet*(15120. + C*h*(840. + C*h*(-210. + C*h*(-15. + 1.875*C*h)))),
      idet*(15120. + C*h*(-3230.83 + C*h*(111.847 + C*h*(22.8034 - 1.11962*C*h)))),
      idet*(15120. + C*h*(-6010.72 + C*h*(961.4 + C*h*(-68.027 + 0.644576*C*h))))
    };
  };

  auto time = 0.0;
  solution.push_back(pair_type{time,y0});

  const auto scale = 1;
  const auto dt = h/scale;

  const auto C = 15.0;
  for (int i = 0; i < nstep*scale; i++)
  {
    const auto& u0 = solution.back();

    const auto solvec = computeSolvec(C,dt);
    auto dy = 0.0;
    for (int i = 0; i < (int)solvec.size(); i++)
    {
      const auto coeff = solvec[i] * get<1>(u0) * dt;
      dy += weights[i]*coeff;
    }
    dy *= -C;

    solution.emplace_back(get<0>(u0)+dt,get<1>(u0)+dy);
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

  const auto h = 1.0/16;

  const auto tmax = 1.0;

  const auto y0 = 1.0;
  const auto solution = solve5(h,tmax, y0);

  for (const auto &x : solution)
  {
//    fprintf(cout, "% %\n", get<0>(x), get<1>(x));
    printf(" %15.15g %15.15g \n", get<0>(x), get<1>(x));
  }


  return 0;

}
