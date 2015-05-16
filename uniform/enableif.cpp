#include <iostream>
#include <utility>
#include <type_traits>

template<typename real_t, typename = void>
struct dtype
{
  static constexpr auto str = "non_floating_point";
};
template<typename real_t>
struct dtype<real_t, typename std::enable_if<std::is_same<real_t,float>::value>::type>
{
  static constexpr auto str = "float";
};
template<typename real_t>
struct dtype<real_t, typename std::enable_if<std::is_same<real_t,double>::value>::type>
{
  static constexpr auto str = "double";
};


int main()
{

  std::cout << dtype<float>::str << std::endl;
  std::cout << dtype<size_t>::str << std::endl;
  std::cout << dtype<double>::str << std::endl;
  std::cout << dtype<int>::str << std::endl;

  return 0;
}
