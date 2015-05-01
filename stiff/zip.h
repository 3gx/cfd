#pragma once
#include<tuple>
#include<iterator>
#include<utility>

/*******************************
 * zipper: adapted from http://codereview.stackexchange.com/questions/30846/zip-like-functionality-with-c11s-range-based-for-loop
 ******************************/

template< bool B, class T = void >
using enableIf = typename std::enable_if<B,T>::type;
using std::index_sequence;
using std::make_index_sequence;

/**************************************
 * increment everu element of a tuple *
 **************************************/
template<size_t... I, typename... Tp>
static inline void increment_tuple_impl(std::tuple<Tp...>& t, index_sequence<I...>)
{
  std::make_tuple( std::get<I>(t)++...  );
}
template<typename... Ts>
static inline void increment_tuple(std::tuple<Ts...>& ts)
{
  increment_tuple_impl(ts,make_index_sequence<sizeof...(Ts)>());
}

/*****************************
 * check equality of a tuple *
 *****************************/
template<typename... Ts>
static inline bool not_equal_tuple_impl(const std::tuple<Ts...>& t1, const std::tuple<Ts...>& t2, index_sequence<>) 
{
  return true;
}
template<size_t I, size_t... Is, typename... Ts>
static inline bool not_equal_tuple_impl(const std::tuple<Ts...>& t1, const std::tuple<Ts...>& t2, index_sequence<I, Is...>) 
{
  return (std::get<I>(t1) != std::get<I>(t2)) && not_equal_tuple_impl(t1,t2,index_sequence<Is...>());
}
template<typename... Ts>
static inline bool not_equal_tuple(const std::tuple<Ts...>& t1, const std::tuple<Ts...> &t2)
{
  return not_equal_tuple_impl(t1,t2, make_index_sequence<sizeof...(Ts)>());
}


/**************************** 
// dereference every element of a tuple (applying operator* to each element, and returning the tuple)
 ****************************/
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

template<typename... Ts >
class zipper
{
  public:

    class iterator : std::iterator<std::forward_iterator_tag, std::tuple<typename Ts::value_type...> >
  {
    protected:
      std::tuple<typename Ts::iterator...> current;
    public:

      explicit iterator(typename Ts::iterator... ts ) : 
        current(ts...) {};

      iterator( const iterator& rhs ) :  current(rhs.current) {};

      iterator& operator++() {
        increment_tuple(current);
        return *this;
      }

      iterator operator++(int) {
        auto a = *this;
        increment_tuple(current);
        return a;
      }

      bool operator!=( const iterator& rhs ) {
        return not_equal_tuple(current, rhs.current);
      }

      typename iterator::value_type operator*() {
        return dereference_tuple(current);
      }
  };


    explicit zipper(Ts&... ts):
      begin_(ts.begin()...), 
      end_( ts.end()...) {};

    zipper(const zipper<Ts...>& a) :
      begin_(  a.begin_ ), 
      end_( a.end_ ) {};

    template<typename... Us>
      zipper<Us...>& operator=( zipper<Us...>& rhs) {
        begin_ = rhs.begin_;
        end_ = rhs.end_;
        return *this;
      }

    zipper<Ts...>::iterator& begin() {
      return begin_;
    }

    zipper<Ts...>::iterator& end() {
      return end_;
    }

    zipper<Ts...>::iterator begin_;
    zipper<Ts...>::iterator end_;
};



//from cppreference.com: 
template <class T>
struct special_decay
{
  using type = typename std::decay<T>::type;
};

//allows the use of references:
template <class T>
struct special_decay<std::reference_wrapper<T>>
{
  using type = T&;
};

template <class T>
using special_decay_t = typename special_decay<T>::type;

//allows template type deduction for zipper:
  template <class... Types>
zipper<special_decay_t<Types>...> zip(Types&&... args)
{
  return zipper<special_decay_t<Types>...>(std::forward<Types>(args)...);
}
