#pragma once

#include <memory>

#include "rti/utils.hpp"

namespace rti {
  template<typename T>
  class i_origin { // Interface
  public:
    virtual ~i_origin() {}
    // Specifies an origin as a subset of R^3.
    // Gives a way to obtain a location from the origin (e.g. in a random manner)
    // Clone function, copies the object
    virtual std::unique_ptr<i_origin<T> > clone() const = 0;
    virtual rti::triple<T> get() = 0;
  };
} // namespace rti
