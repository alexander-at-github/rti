#pragma once

#include "rti/utils.hpp"

namespace rti {
  // Interface
  template<typename Ty>
  class i_abs_geometry {
  public:
    // virtual destructor
    virtual ~i_abs_geometry() {}
    virtual void print(std::ostream& pOs) const = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual rti::triple<Ty> get_normal(unsigned int primID) const = 0;
    virtual rti::triple<Ty> get_new_origin(unsigned int primID) const = 0;
  };
} // namespace rti
