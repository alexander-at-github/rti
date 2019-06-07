#pragma once

#include <memory>

#include "rti/types.hpp"

namespace rti {
  template<typename T>
  class i_direction { // Interface
    // Specifies a direction (independent of an origin)
    // Gives a way to obtain a concret direction (e.g. in a random manner)
  public:
    virtual std::unique_ptr<i_direction<T> > clone() const = 0;
    virtual rti::triple<T> get() = 0;
  };
} // namespace rti
