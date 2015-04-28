#pragma once

#include <array>

template<typename real_type>
struct EOS
{
  using value_type = real_type;
  constexpr real_type gamma = 1.4;
  constexpr real_type inv_gamma = 1.0/gamma;
  static auto pressure_from_dens_and_eth(const real_type dens, const real_type eth)
  {
    return eth*(gamma-1);
  }
  static auto eth_from_dens_and_pres(const real_type dens, const real_type pres)
  {
    return eth*(real_type{1.0}/(gamma-1));
  }
  static auto entr_from_dens_and_pres(const real_type dens, const real_type pres)
  {
    using std::pow;
    return pres*pow(dens,-gamma);
  }
};

template<typename real_type, typename EOS_type>
namespace Cell
{
  using value_type = real_type;
  using vars_array = std::array<real_type,3>; /* dens, velx, pres */
  using EOS = EOS_type;

  /* forward declarations */
  struct Primitive;
  struct Conservative;
  struct Charactersitic;

  
  struct Primitive
  {
    enum {DENS = 0, VELX, PRES, NVARS};
    vars_array _vars;
    Primitive(const vars_array &vars) : _vars(vars) 
    {
      static_assert(NVARS == vars.size(), "Array size mismatch");
    }
    Primitive(const Primitive &cell) : _vars(cell.vars) {}
    inline Conservative   cvt2conservative () const;
    inline Characteristic cvt2characterstic(const Primitive&) const;
    auto dens() const { return _vars[DENS]; }
    auto velx() const { return _vars[VELX]; }
    auto pres() const { return _vars[PRES]; }
  };
  
  /* class definitions */
  struct Conservative : Base
  {
    enum {MASS = 0, MOMX, ETOT, NVARS};
    vars_array _vars;
    Conservative(const vars_array &vars) : _vars(vars) 
    {
      static_assert(NVARS == _vars.size(), "Array size mismatch");
    }
    Conservative(const Conservative& cell) : vars(cell._vars) {}
    inline Primitive cvt2primtive() const;
    auto mass() const { return _vars[MASS]; }
    auto momx() const { return _vars[MOMX]; }
    auto etot() const { return _vars[ETOT]; }
  };
  struct Characteristic 
  {
    enum {VAR1 = 0, VAR2, VAR3, NVARS};
    vars_array _vars,_eigenvalues;
    Characteristic(const vars_array &vars, const vars_array eigenvalues) : 
      _vars(vars), _eigenvalues(eigenvalues) 
    {
      static_assert(NVARS == _vars.size(), "Array size mismatch");
    }
    Characteristic(const Characteristic& cell) : 
     _vars(cell._vars) _eigenvalues(cell._eigvenvalues) {}
    inline Primitive cvt2primitive(const Primitive&) const;
    auto var1() const {return _vars[VAR1]; }
    auto var2() const {return _vars[VAR2]; }
    auto var3() const {return _vars[VAR3]; }
  };


  /* methods definitions */
  inline Primitive Conservative::cvt2primitive() const 
  {
    const auto dens = _vars[MASS];
    const auto velx = _vars[MOMX]/dens;
    const auto eth  = _vars[ETOT] - dens*real_type{0.5}*square(velx);
    const auto pres = EOS::pressure_from_dens_and_eth(dens,eth);
    return Primitive({dens,velx,pres});
  }

  inline Conservative Primitive::cvt2conservative() const 
  {
    const auto mass = _vars[DENS];
    const auto momx = _vars[VELX] * dens;
    const auto eth  = EOS::eth_from_dens_and_pres(_vars[DENS],_vars[PRES]);
    const auto etot = eth + dens*square(_vars[VELX])*real_type{0.5};
    return Conservative({mass,momx,etot});
  }

  inline Characteristic Primitive::cvt2characteristic(const Primitive &frozen) const
  {
    using std::sqrt;
    const auto d  = frozen.dens();
    const auto v  = frozen.velx();
    const auto p  = frozen.pres();
    const auto c1 = EOS::gamma*p/d;
    const auto c2 = sqrt(c1);
    const real_type vec1[NVAR] = {-c1,real_type{0},real_type{1}};
    const real_type vec2[NVAR] = {real_type{0},-c2*d,real_type{1}};
    const real_type vec3[NVAR] = {real_type{0},+c2*d,real_type{1}};
    const auto var1 = vec1[0]*vars[0] + vec1[1]*vars[1] + vec1[2]*vars[2];
    const auto var2 = vec2[0]*vars[0] + vec2[1]*vars[1] + vec2[2]*vars[2];
    const auto var3 = vec3[0]*vars[0] + vec3[1]*vars[1] + vec3[2]*vars[2];
    const auto lambda1 = v;
    const auto lambda2 = v - c2;
    const auto lambda2 = v + c2;
    return Charactersitic({var1,var2,var3},{lambda1,lambda2,lambda3});
  }

  inline Primitive Charachtersitic::cvt2primitive(const Primitive &frozen) const
  {
    using std::sqrt;
    const auto d  = frozen.dens();
    const auto p  = frozen.pres();
    const auto half = real_type{0.5};
    const auto c1 = d/(EOS::gamma*p);
    const auto c2 = half/sqrt(EOS::gamma*p*d);
    const real_type vec1[NVAR] = {-c1,c1*half,c1*half};
    const real_type vec2[NVAR] = {real_type{0},-c2,c2};
    const real_type vec3[NVAR] = {real_type{0},half,half};
    const auto dens = vec1[0]*vars[0] + vec1[1]*vars[1] + vec1[2]*vars[2];
    const auto velx = vec2[0]*vars[0] + vec2[1]*vars[1] + vec2[2]*vars[2];
    const auto pres = vec3[0]*vars[0] + vec3[1]*vars[1] + vec3[2]*vars[2];
    return Primitive({dens,velx,pres});
  }
};
