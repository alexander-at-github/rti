#pragma once

#include <memory>

#include "rti/util/utils.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class i_source { // Interface
  public:
    virtual ~i_source() {}
    // Takes a RTCRay class and sets its members, e.g., origin and direction
    virtual void fill_ray(RTCRay&, rti::rng::i_rng&, rti::rng::i_rng::i_state&) = 0;
    virtual void optional_consider(rti::util::triple<Ty> SourceXyz, double relativeerror) = 0;
    virtual void optional_update_adaptive_sampling() = 0;
  };
}}
