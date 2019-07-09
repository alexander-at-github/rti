#pragma once

#include <memory>

#include "rti/utils.hpp"

namespace rti {
  template<typename T>
  class i_direction { // Interface
    // Specifies a direction (independent of an origin)
    // Gives a way to obtain a concret direction (e.g. in a random manner)
  public:
    virtual ~i_direction() {}
    virtual rti::triple<T> get(rti::i_rng&, rti::i_rng::i_state&) const = 0;
  };
} // namespace rti
