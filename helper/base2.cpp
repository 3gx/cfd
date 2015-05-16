#include <iostream>

constexpr int Log2(int n, int p = 0) 
{
  return (n <= 1) ? p : Log2(n / 2, p + 1);
}
template<int N>
constexpr int Log2()
{
  static_assert((1 << Log2(N,0)) == N, "N is not a power of 2!");
  return Log2(N,0);
}

int main()
{
  std::cout << Log2<32>() << std::endl;
};


