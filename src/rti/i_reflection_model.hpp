#pragma once

namespace rti {
  class i_reflection_model {
    // Pure Virtual Class
    virtual ~i_reflection_model() {};
    virtual bool use(const RTCRayHit& pRayhit) = 0;
  };
} // namespace rti
