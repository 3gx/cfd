#include <iostream>


template<typename parms_t, typename vector_t>
void advection_flux(const parms_t &parms, const vector_t &f, vector_t &flux)
{
  using real_t = vector_t::value_type;
  const auto n = state.size();
  for (int i = 0; i < n; i++)
  {
    const auto b1 = 
      13.0/12.0 * square(f[i-2] - 2.0*f[i-1] +     f[i]) +
       1.0/ 4.0 * square(f[i-2] - 4.0*f[i-1] + 3.0*f[i]);

    const auto b2 =
      (13.0/12.0) * square(f[i-1] - 2.0*f[i] + f[i+1]) +
      ( 1.0/ 4.0) * square(f[i-1]            + f[i+1]);

    const auto b3 = 
      (13.0/12.0) * square(    f[i] - 2.0*f[i+1] + f[i+2]) +
      ( 1.0/ 4.0) * square(3.0*f[i] - 4.0*f[i+1] + f[i+2]);

    const auto g1 = real_t{0.1};
    const auto g2 = real_t{0.6};

    auto wh1 = g1/std::pow(WENO_EPS + b1, WENO_POW);
    auto wh2 = g2/std::pow(WENO_EPS + b1, WENO_POW);
    auto wh3 = g3/std::pow(WENO_EPS + b1, WENO_POW);

    auto wh = 1.0/(wh1 + wh2 + wh3);
    wh1 *= wh;
    wh2 *= wh;
    wh3 *= wh;

    auto fh1 = ( 1.0/3.0)*f[i-2] - (7.0/6.0)*f[i-1] + (11.0/6.0)*f[i  ];
    auto fh2 = (-1.0/6.0)*f[i-1] + (5.0/6.0)*f[i  ] + ( 1.0/3.0)*f[i+1];
    auto fh3 = ( 1.0/3.0)*f[i  ] + (5.0/6.0)*f[i+1] + (-1.0/6.0)*f[i+1];

    auto flux = wh1*fh1 + wh2*fh2 + wh3*fh3;


  }
};

template<typename vector_t>
void diffusion_flux(const vector_t &state, vector_t &flux)
{
}

int main(int argc, char *argv[])
{
  return 0;
}
