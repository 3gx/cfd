#pragma once

#include "common.h"

template<typename real_type, typename EOS_type>
struct WenoFD5
{
  using value_type = real_type;
  using EOS = EOS_type;

  using Cell = CellT<real_type,EOS>;

  std::vector<Cell::Conservative> cons_vec;

};


