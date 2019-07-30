#pragma once

#include "rti/ray/i_direction.hpp"
#include "rti/ray/i_origin.hpp"
#include "rti/ray/i_source.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class source : public rti::ray::i_source {
    // Combines an origin with a direction
    // Not thread safe. See direction classes for reasons.
  public:

    source(rti::ray::i_origin<Ty>& pOrigin, rti::ray::i_direction<Ty>& pDirection) :
      mOrigin(pOrigin),
      mDirection(pDirection) {}

    void fill_ray(RTCRay& pRay, rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final {
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
    rti::ray::i_origin<Ty>& mOrigin;
    rti::ray::i_direction<Ty>& mDirection;
  };
}} // namespace
