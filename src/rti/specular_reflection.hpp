#pragma once

#include "rti/i_reflection_model.hpp"

namespace rti {
  class specular_reflection : public i_reflection_model {
  public:
    bool use(const RTCRayHit& pRayhit) const override final {
      assert(false && "Not implemented");
      return false;
    }
    // TODO: create modular random number abstraction
  };
} // namespace rti
