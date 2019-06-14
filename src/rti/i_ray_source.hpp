#pragma once

#include <memory>

namespace rti {
  class i_ray_source { // Interface
  public:
    virtual ~i_ray_source() {}
    virtual std::unique_ptr<i_ray_source> clone() const = 0;
    // Takes a RTCRay class and sets its members, e.g., origin and direction
    virtual void fill_ray(RTCRay&) = 0;
  };
} // namespace rti
