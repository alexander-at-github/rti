#pragma once

#include <memory>

#include "rti/rng/i_rng.hpp"
#include "rti/util/utils.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class i_origin { // Interface
  public:
    virtual ~i_origin() {}
    // Specifies an origin as a subset of R^3.
    // Gives a way to obtain a location from the origin (e.g. in a random manner)
    virtual rti::util::triple<Ty> get(rti::rng::i_rng&, rti::rng::i_rng::i_state&) = 0;
  };
}} // namespace
