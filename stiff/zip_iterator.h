#pragma once
#include<tuple>
#include<iterator>
#include<utility>

/*******************************
 * ZipIterator: 
 * adapted from http://codereview.stackexchange.com/questions/30846/zip-like-functionality-with-c11s-range-based-for-loop
 ******************************/

template< bool B, class T = void >
using _enableIf = typename std::enable_if<B,T>::type;
using std::index_sequence;
using std::make_index_sequence;

/**************************************
 * increment every element of a tuple *
 **************************************/

  template<typename Tuple, size_t... I>
static inline void increment_tuple_impl(Tuple& t, index_sequence<I...>)
{
  std::make_tuple( ++std::get<I>(t)...  );
}
template<typename... Ts>
static inline void increment_tuple(std::tuple<Ts...>& ts)
{
  increment_tuple_impl(ts,make_index_sequence<sizeof...(Ts)>());
}

/********************************
 * check equality of two tuples *
 ********************************/
template<typename Tuple>
static inline bool not_equal_tuple_impl(const Tuple& t1, const Tuple& t2, index_sequence<>)
{
  return true;
}
template<typename Tuple, size_t I, size_t... Is>
static inline bool not_equal_tuple_impl(const Tuple& t1, const Tuple& t2, index_sequence<I, Is...>)
{
  return (std::get<I>(t1) != std::get<I>(t2)) && not_equal_tuple_impl(t1,t2,index_sequence<Is...>());
}
template<typename... Ts>
static inline bool not_equal_tuple(const std::tuple<Ts...>& t1, const std::tuple<Ts...> &t2)
{
  return not_equal_tuple_impl(t1,t2, make_index_sequence<sizeof...(Ts)>());
}


/****************************************
 * dereference every element of a tuple *
 ****************************************/
template <typename Tuple, size_t... I>
static inline auto dereference_tuple_impl(const Tuple& t, index_sequence<I...>) ->
decltype(std::tie(*std::get<I>(t)...))
{
  return std::tie(*std::get<I>(t)...);
}
template<typename... Ts>
static inline auto dereference_tuple(std::tuple<Ts...>& ts) -> 
decltype( dereference_tuple_impl( std::tuple<Ts...>(), make_index_sequence<sizeof...(Ts)>()))
{
  return dereference_tuple_impl( ts, make_index_sequence<sizeof...(Ts)>());
}

template<typename... Types>
class ZipIterator
{
  public:

    template<typename T>
      using remove_ref = typename std::remove_reference<T>::type;

    template<bool Condition, typename Then, typename Else>
      struct _if { using type = Then; };
    template<typename Then, typename Else>
      struct _if<false,Then,Else> { using type = Else; };

    template<typename Container>
      using iterator_type = typename _if<
            std::is_const<remove_ref<Container>>::value, 
            typename remove_ref<Container>::const_iterator,
            typename remove_ref<Container>::iterator > :: type;

    template<typename Type>
      using value_type = typename _if<
            std::is_const<remove_ref<Type>>::value, 
            typename std::add_const<typename remove_ref<Type>::value_type>::type,
            typename remove_ref<Type>::value_type > :: type;

    class Iterator : std::iterator<std::forward_iterator_tag, std::tuple<value_type<Types>&...>>
  {
    protected:
      std::tuple<iterator_type<Types>...> current;
    public:

      explicit Iterator(iterator_type<Types>... ts ) :
        current(ts...) {};

      Iterator( const Iterator& rhs ) :  current(rhs.current) {};

      Iterator& operator++() {
        increment_tuple(current);
        return *this;
      }

      auto operator++(int) {
        auto a = *this;
        increment_tuple(current);
        return a;
      }

      auto operator!=( const Iterator& rhs ) {
        return not_equal_tuple(current, rhs.current);
      }

      auto operator*()
      {
        return dereference_tuple(current);
      }
  };


    explicit ZipIterator(Types&&... ts):
      _begin(ts.begin()...), 
      _end( ts.end()...)  {}

    Iterator& begin() { return _begin; }
    Iterator& end  () { return _end; }

    Iterator _begin;
    Iterator _end;
};

template <class... Types>
auto make_zip_iterator(Types&&... ts)
{
  return ZipIterator<Types...>(std::forward<Types>(ts)...);
}
