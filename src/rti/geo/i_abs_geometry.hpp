#pragma once

#include "rti/util/utils.hpp"

namespace rti { namespace geo {
  // Interface
  template<typename Ty>
  class i_abs_geometry {
  public:
    // virtual destructor
    virtual ~i_abs_geometry() {}
    virtual void print(std::ostream& pOs) const = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual rti::util::triple<Ty> get_normal(unsigned int primID) const = 0;
    virtual rti::util::triple<Ty> get_new_origin(RTCRayHit& pRayhit, unsigned int primID) const = 0;
  };
}} // namespace
