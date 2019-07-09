#pragma once

#include <memory>

namespace rti {
  class i_ray_source { // Interface
  public:
    virtual ~i_ray_source() {}
    // Takes a RTCRay class and sets its members, e.g., origin and direction
    virtual void fill_ray(RTCRay&, rti::i_rng&, rti::i_rng::i_state&) const = 0;
  };
} // namespace rti
