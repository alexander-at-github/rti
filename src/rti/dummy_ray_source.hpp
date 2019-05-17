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
      ray.dir_x = this->get_num();
      ray.dir_y = this->get_num();
      ray.dir_z = this->get_num();
      // start of ray
      ray.tnear = 0;
      // Maximum length of ray
      ray.tfar = std::numeric_limits<float>::max();
      BOOST_LOG_SEV(rti::mRLogger, blt::debug) << "WARNING: This is a dummy ray source.";
      return ray;
    }
  private:
    float get_num() {
      ///////////////////////////////////////////////////////////////
      // This is an entirely uninformed pseudo random float generator
      // in the interval [0,1].
      ///////////////////////////////////////////////////////////////
      static int xx = 11111111;
      xx ^= xx << 2;
      xx += 1;
      return (static_cast<float>(xx) / std::numeric_limits<int>::max());
    }
  };
}
