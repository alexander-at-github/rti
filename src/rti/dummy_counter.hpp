#pragma once

#include "rti/i_hit_counter.hpp"

namespace rti {
  class dummy_counter : public i_hit_counter {
  public:
    void use(const RTCRayHit& pRayhit) override final {
      return;
    }
  };
} // namespace rti
