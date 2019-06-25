#pragma once

#include "rti/i_reflection_model.hpp"

namespace rti {
  class specular_reflection : public i_reflection_model {
  public:
    // Decides whether or not to reflect. If a reflection should happen, it sets the origin and direction
    // in pRayhit and returns true. If no reflection should happen, then it does not change pRayhit and
    // returns false.
    bool use(const RTCRayHit& pRayhit) override final {
      assert(false && "Not implemented");
      return false;
    }
  };
} // namespace rti
