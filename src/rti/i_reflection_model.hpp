#pragma once

namespace rti {
  class i_reflection_model {
  public:
    // Pure Virtual Class
    virtual ~i_reflection_model() {};
    // Decides whether or not to reflect. If a reflection should happen, it sets the origin and direction
    // in pRayhit and returns true. If no reflection should happen, then it does not change pRayhit and
    // returns false.
    virtual bool use(const RTCRayHit& pRayhit) const = 0;
  };
} // namespace rti
