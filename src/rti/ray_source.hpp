#pragma once

#include "rti/i_direction.hpp"
#include "rti/i_origin.hpp"
#include "rti/i_ray_source.hpp"

namespace rti {
  template<typename T>
  class ray_source : public i_ray_source {
    // Combines an origin with a direction
    // Not thread safe. See direction classes for reasons.
  public:

    ray_source(rti::i_origin<T>& pOrigin, rti::i_direction<T>& pDirection) :
      mOrigin(pOrigin),
      mDirection(pDirection) {}

    void fill_ray(RTCRay& pRay, rti::i_rng& pRng, rti::i_rng::i_state& pRngState) const override final {
      auto origin = mOrigin.get(pRng, pRngState);
      pRay.org_x = (float) origin[0];
      pRay.org_y = (float) origin[1];
      pRay.org_z = (float) origin[2];
      auto direction = mDirection.get(pRng, pRngState);
      pRay.dir_x = (float) direction[0];
      pRay.dir_y = (float) direction[1];
      pRay.dir_z = (float) direction[2];
    }

  private:
    i_origin<T>& mOrigin;
    i_direction<T>& mDirection;
  };
} // namespace rti
