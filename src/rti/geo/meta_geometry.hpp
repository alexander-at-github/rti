#pragma once

#include "../util/utils.hpp"

namespace rti { namespace geo {
  // Interface
  template<typename Ty>
  class meta_geometry {
  public:
    // virtual destructor
    virtual ~meta_geometry() {}
    virtual void print(std::ostream& pOs) = 0;
    virtual RTCDevice& get_rtc_device() = 0;
    virtual RTCGeometry& get_rtc_geometry() = 0;
    virtual rti::util::triple<Ty> get_normal(unsigned int primID) = 0;
    virtual rti::util::triple<Ty>& get_normal_ref(unsigned int primID) = 0;
    virtual rti::util::triple<Ty> get_new_origin(RTCRay& pRay, unsigned int primID) 
    {
      auto xx = pRay.org_x + pRay.dir_x * pRay.tfar;
      auto yy = pRay.org_y + pRay.dir_y * pRay.tfar;
      auto zz = pRay.org_z + pRay.dir_z * pRay.tfar;
      return {(Ty) xx, (Ty) yy, (Ty) zz};
    }
  };
}}
