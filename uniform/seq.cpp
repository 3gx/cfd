#include <iostream>
#include <utility>

template< bool B, class T = void >
using enableIf = typename std::enable_if<B,T>::type;

template<size_t... I>
using indexSeq = std::index_sequence<I...>;

template<size_t N>
using makeIndexSeq = std::make_index_sequence<N>;

template<size_t N, size_t... S>
struct makeZeroSeqImpl : makeZeroSeqImpl<N-1,0,S...> {};
template<size_t... S>
struct makeZeroSeqImpl<0,S...>
{
  using type = indexSeq<S...>;
};

template<size_t N>
using makeZeroSeq = typename makeZeroSeqImpl<N>::type;


template<size_t... I>
void foo(indexSeq<I...>)
{
  std::cout << __PRETTY_FUNCTION__ << std::endl;
}

int main()
{
  foo(makeIndexSeq<3>());
  foo(makeZeroSeq<5>());
}
