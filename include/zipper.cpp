#include <iostream>
#include <vector>
#include <tuple>
#include "zip_iterator.h"

template<size_t... I>
void foo(index_sequence<I...>)
{
  std::cout << __PRETTY_FUNCTION__<< std::endl;
}

#if 0
  template<size_t... I>
void goo(ct_integers_list<I...>)
{
  std::cout << __PRETTY_FUNCTION__<< std::endl;
}
#endif

template<typename... Ts>
void bar(Ts... ts)
{
  foo(make_index_sequence<sizeof...(Ts)>());
  //goo(typename ct_iota_1<sizeof...(Ts)>::type());
}

template<typename... TYPE>
struct TD
{
  static void foo()
  {
    std::cout << __PRETTY_FUNCTION__<<  std::endl;
  }

};

template<typename... Ts>
struct zips
{
  explicit zips(Ts&& ... ts)
  {

  };
};
template<typename... Ts>
zips<Ts&...> ziptest(Ts&&... ts)
{
  TD<typename Ts::iterator...>::foo();
  return zips<Ts...>(std::forward<Ts>(ts)...);
}

int main()
{
  bar(1.0,2U,'c');
  using pair_type = std::tuple<int,float>;
  using std::get;

  const std::vector<int> ai{1,2,3,4,5};
  std::vector<float> af{1.1,2.2,3.3,4.4,5.5};
  std::vector<pair_type> aif(af.size(),std::make_tuple(-1,-1.0));


#if 0
  for (size_t i = 0; i < aif.size(); i++)
  {
    aif[i] = std::tie(ai[i],af[i]);
  }
#endif
  for (const auto& x : ai)
  {
    std::cout << x << std::endl;
  }
  for (auto x : af)
  {
    std::cout << x << std::endl;
    x=3;
  }

  std::cout << " ----- \n";

#if 0
  for (auto val : zip(ai,af))
  {
    std::cout << std::get<0>(val) << " - " << std::get<1>(val) << std::endl;
  }
#endif
#if 1
  auto zipit = make_zip_iterator(aif,ai,af);
  auto zipit1 = make_zip_iterator(aif,ai,af);
  zipit1 = zipit;
  for (auto val : zipit1) //make_zip_iterator(aif,ai,af))
  {
    std::get<0>(val) = std::make_tuple(std::get<1>(val), std::get<2>(val));
  }
#endif
  std::cout << " ----- \n";
  for (const auto& val : aif)
  {
    std::cout << std::get<0>(val) << " - " << std::get<1>(val) << std::endl;
  }

}
