#pragma once

#include <memory>

namespace rti { namespace ray {
  class i_source { // Interface
  public:
    virtual ~i_source() {}
    // Takes a RTCRay class and sets its members, e.g., origin and direction
    virtual void fill_ray(RTCRay&, rti::rng::i_rng&, rti::rng::i_rng::i_state&) const = 0;
  };
}} // namespace
