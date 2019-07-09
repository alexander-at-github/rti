#pragma once

#include <memory>

#include "rti/i_rng.hpp"
#include "rti/utils.hpp"

namespace rti {
  template<typename T>
  class i_origin { // Interface
  public:
    virtual ~i_origin() {}
    // Specifies an origin as a subset of R^3.
    // Gives a way to obtain a location from the origin (e.g. in a random manner)
    virtual rti::triple<T> get(rti::i_rng&, rti::i_rng::i_state&) const = 0;
  };
} // namespace rti
