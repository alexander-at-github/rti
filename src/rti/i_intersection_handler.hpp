#pragma once

namespace rti {
  class i_intersection_handler {
  public:
    // Pure Virtual Class
    virtual ~i_intersection_handler() {};
    virtual void use(const RTCRayHit& pRayhit) = 0;
  };
} // namespace rti
