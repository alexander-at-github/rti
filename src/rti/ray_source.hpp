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

    ray_source(
        std::unique_ptr<rti::i_origin<T> > pOrigin,
        std::unique_ptr<rti::i_direction<T> > pDirection) :
      mOrigin(std::move(pOrigin)),
      mDirection(std::move(pDirection)) {}

    // Note: With unique_ptr one cannot return covariant types
    std::unique_ptr<i_ray_source> clone() const override final {
      RLOG_TRACE << "[ray_source.clone()]" << std::endl;
      return std::make_unique<ray_source<T> >(mOrigin->clone(), mDirection->clone());
    }

    void fill_ray(RTCRay& pRay) override final {
      auto origin = mOrigin->get();
      pRay.org_x = (float) origin.frst;
      pRay.org_y = (float) origin.scnd;
      pRay.org_z = (float) origin.thrd;
      auto direction = mDirection->get();
      pRay.dir_x = (float) direction.frst;
      pRay.dir_y = (float) direction.scnd;
      pRay.dir_z = (float) direction.thrd;
    }

  private:
    std::unique_ptr<i_origin<T> > mOrigin;
    std::unique_ptr<i_direction<T> > mDirection;
  };
} // namespace rti
