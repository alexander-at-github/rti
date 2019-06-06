#pragma once

//#include<embree3/rtcore.h>

namespace rti {
  class i_ray_source {
  public:
    // Interface
    virtual~i_ray_source() {}
    virtual RTCRay get_ray() = 0;
    // Takes a RTCRay class and sets its members, e.g., origin and direction
    virtual void set_ray(RTCRay&) = 0;
  };
} // namespace rti
