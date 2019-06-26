#pragma once

#include "rti/i_reflection_model.hpp"

namespace rti {
  class lambertian_reflection : public i_reflection_model {
  public:
    bool use(RTCRayHit& pRayhit, const i_geometry& pGeometry, i_hit_counter& pHitcounter) const override final {
      // TODO: Get random number and decide whether or not to reflect.
      // TODO: get surface normal at intersection
      // TODO: Compute lambertian reflection with respect to surface normal

      // TODO: FIX
      pHitcounter.use(pRayhit);
      return false;
   }
  };
} // namespace rti
