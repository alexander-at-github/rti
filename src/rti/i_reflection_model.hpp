#pragma once

namespace rti {
  class i_reflection_model {
  public:
    // Pure Virtual Class
    virtual ~i_reflection_model() {}
    // Decides whether or not to reflect. If a reflection should happen, it sets
    // the origin and direction in the RTCRayHit object and returns true. If no
    // reflection should happen, then it does not change pRayhit and returns
    // false.
    virtual bool use(RTCRayHit&, rti::i_rng&, rti::i_rng::i_state&, const i_geometry&, i_hit_counter&) const = 0;
  };
} // namespace rti
