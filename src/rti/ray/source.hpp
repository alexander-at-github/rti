#pragma once

#include <x86intrin.h> // vector instruction instrinsics

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

      // "Avoid store-to-load forwarding issues with single rays
      //
      // We recommend to use a single SSE store to set up the org and tnear components,
      // and a single SSE store to set up the dir and time components of a single ray (RTCRay type).
      // Storing these values using scalar stores causes a store-to-load forwarding penalty because
      // Embree is reading these components using SSE loads later on." Source: https://www.embree.org/api.html

      auto tnear = 1e-4f; // float

      auto orgn = mOrigin.get(pRng, pRngState);
      // pRay.org_x = (float) orgn[0];
      // pRay.org_y = (float) orgn[1];
      // pRay.org_z = (float) orgn[2];
      // pRay.tnear = tnear;

      // float vara[4] = {(float) orgn[0], (float) orgn[1], (float) orgn[2], tnear};
      // reinterpret_cast<__m128&>(pRay) = _mm_load_ps(vara);

      // the following instruction would have the same result
      // the intrinsic _mm_set_ps turns the ordering of the input around.
      reinterpret_cast<__m128&>(pRay) = _mm_set_ps(tnear, (float) orgn[2], (float) orgn[1], (float) orgn[0]);

      auto time = 0.0f; // float

      auto dir = mDirection.get(pRng, pRngState);
      // pRay.dir_x = (float) dir[0];
      // pRay.dir_y = (float) dir[1];
      // pRay.dir_z = (float) dir[2];
      // pRay.time = time;

      // float varb[4] = {(float) dir[0], (float) dir[1], (float) dir[2], time};
      // reinterpret_cast<__m128&>(pRay.dir_x) = _mm_load_ps(varb);
      reinterpret_cast<__m128&>(pRay.dir_x) = _mm_set_ps(time, (float) dir[2], (float) dir[1], (float) dir[0]);
    }

    rti::ray::i_origin<Ty>&
    get_origin()
    {
      return mOrigin;
    }
    rti::ray::i_direction<Ty>&
    get_direction()
    {
      return mDirection;
    }


  private:
    rti::ray::i_origin<Ty>& mOrigin;
    rti::ray::i_direction<Ty>& mDirection;
  };
}} // namespace
