#pragma once

#include "i_direction.hpp"
#include "../rng/i_rng.hpp"
#include "../util/utils.hpp"

namespace rti { namespace ray {
  template<typename numeric_type>
  // cosine direction in opposite direction of z achsis
  class cosine_direction_z : public ray::i_direction<numeric_type> {

  public:
    util::triple<numeric_type>
    get(rng::i_rng& pRng, rng::i_rng::i_state& pRngState) const override final
    {
      auto r1 = ((numeric_type) pRng.get(pRngState)) / ((numeric_type) pRng.max());
      auto r2 = ((numeric_type) pRng.get(pRngState)) / ((numeric_type) pRng.max());
      auto zz = - sqrtf(r2);
      auto xx = cosf(two_pi * r1) * sqrtf(1 - r2);
      auto yy = sinf(two_pi * r1) * sqrtf(1 - r2);
      auto result = util::triple<numeric_type> {xx, yy, zz};
      util::normalize(result);
      return result;
    }

  private:
    constexpr static numeric_type two_pi = util::pi() * 2;
  };
}}
