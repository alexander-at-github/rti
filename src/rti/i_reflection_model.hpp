#pragma once

namespace rti {
  class i_reflection_model {
  public:
    // Pure Virtual Class
    virtual ~i_reflection_model() {};
    // Decides whether or not to reflect. If a reflection should happen, it sets
    // the origin and direction in the RTCRayHit object and returns true. If no
    // reflection should happen, then it does not change pRayhit and returns
    // false.
    virtual bool use(RTCRayHit& pR, const i_geometry& pG, i_hit_counter& pH) const = 0;
  };
} // namespace rti
