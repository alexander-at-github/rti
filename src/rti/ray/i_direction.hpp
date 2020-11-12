#pragma once

#include <memory>

#include "../util/utils.hpp"
#include "../rng/i_rng.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class i_direction { // Interface
    // Specifies a direction (independent of an origin)
    // Gives a way to obtain a concret direction (e.g. in a random manner)
  public:
    virtual ~i_direction() {}
    virtual rti::util::triple<Ty> get(rti::rng::i_rng&, rti::rng::i_rng::i_state&) const = 0;
  };
}} // namespace
