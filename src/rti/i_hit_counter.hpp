#pragma once

namespace rti {
  class i_hit_counter {
  public:
    // Pure Virtual Class
    virtual ~i_hit_counter() {};
    virtual void use(const RTCRayHit& pRayhit) = 0;
  };
} // namespace rti
