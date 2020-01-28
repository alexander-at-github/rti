#pragma once

#include "rti/ray/adaptive_origin.hpp"
#include "rti/ray/i_source.hpp"
#include "rti/ray/source.hpp"

namespace rti { namespace ray {
  template<typename Ty>
  class adaptive_source : public rti::ray::i_source<Ty> {
  public:
    adaptive_source(rti::ray::adaptive_origin<Ty>& pOrigin, rti::ray::i_direction<Ty>& pDirection) :
      origin(pOrigin),
      direction(pDirection) {
      std::cout << "[adaptive_source] this class sets the Embree variable tnear to a constant" << std::endl;
    }

    void fill_ray(RTCRay& pRay, rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {
      auto tnear = 1e-4f;
      auto originval = origin.get(pRng, pRngState);
      reinterpret_cast<__m128&>(pRay) =
        _mm_set_ps(tnear, (float) originval[2], (float) originval[1], (float) originval[0]);

      auto time = 0.0f;
      auto dir = direction.get(pRng, pRngState);
      reinterpret_cast<__m128&>(pRay.dir_x) =
        _mm_set_ps(time, (float) dir[2], (float) dir[1], (float) dir[0]);
    }

    // one should probably pass in the maximum relative error of this ray path
    void optional_consider(rti::util::triple<Ty> xyz, double relativeerror) {
      origin.consider(xyz, relativeerror);
    }

    void optional_update_adaptive_sampling() {
      origin.update_adaptive_sampling_state();
    }

  private:
    rti::ray::adaptive_origin<Ty>& origin;
    rti::ray::i_direction<Ty>& direction;
  };
}}
