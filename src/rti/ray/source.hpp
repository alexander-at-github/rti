#pragma once

#include <x86intrin.h> // vector instruction instrinsics

#include "i_direction.hpp"
#include "i_origin.hpp"
#include "i_source.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class source : public ray::i_source {
    // Combines an origin with a direction
    // Not thread safe. See direction classes for reasons.
  public:

    source(ray::i_origin<Ty>& pOrigin, ray::i_direction<Ty>& pDirection) :
      mOrigin(pOrigin),
      mDirection(pDirection) {
      std::cout << "USING source" << std::endl;
    }

    void fill_ray(RTCRay& pRay, rng::i_rng& pRng,
                  rng::i_rng::i_state& pRngState1, rng::i_rng::i_state& pRngState2,
                  rng::i_rng::i_state& pRngState3, rng::i_rng::i_state& pRngState4
                  ) override final {

      // "Avoid store-to-load forwarding issues with single rays
      //
      // We recommend to use a single SSE store to set up the org and tnear components,
      // and a single SSE store to set up the dir and time components of a single ray (RTCRay type).
      // Storing these values using scalar stores causes a store-to-load forwarding penalty because
      // Embree is reading these components using SSE loads later on." Source: https://www.embree.org/api.html

      auto tnear = 1e-4f; // float

      auto orgn = mOrigin.get(pRng, pRngState1, pRngState2);
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

      auto dir = mDirection.get(pRng, pRngState3, pRngState4);
      // pRay.dir_x = (float) dir[0];
      // pRay.dir_y = (float) dir[1];
      // pRay.dir_z = (float) dir[2];
      // pRay.time = time;

      // float varb[4] = {(float) dir[0], (float) dir[1], (float) dir[2], time};
      // reinterpret_cast<__m128&>(pRay.dir_x) = _mm_load_ps(varb);
      reinterpret_cast<__m128&>(pRay.dir_x) = _mm_set_ps(time, (float) dir[2], (float) dir[1], (float) dir[0]);
    }

  private:
    ray::i_origin<Ty>& mOrigin;
    ray::i_direction<Ty>& mDirection;
  };
}} // namespace
