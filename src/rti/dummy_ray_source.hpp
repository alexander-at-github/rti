#pragma once

#include <embree3/rtcore.h>

#include "rti/i_ray_source.hpp"

namespace rti {
  class dummy_ray_source : public i_ray_source {
  public:
    RTCRay get_ray() override final {
      RTCRay ray;
      // Origin:
      ray.org_x = 1;
      ray.org_y = 0;
      ray.org_z = 0;
      // Direction (may not be 0 0 0)
      ray.dir_x = 1;
      ray.dir_y = 1;
      ray.dir_z = 1;
      // start of ray
      ray.tnear = 0;
      // Maximum length of ray
      ray.tfar = std::numeric_limits<float>::max();
      BOOST_LOG_SEV(rti::mRLogger, blt::warning) << "WARNING: This is a dummy ray source.";
      return ray;
    }
  };
}
