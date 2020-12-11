#pragma once

#include "i_direction.hpp"
#include "../rng/i_rng.hpp"
#include "../util/utils.hpp"

namespace rti { namespace ray {
  template<typename numeric_type>
  // power cosine direction in opposite direction of z achsis
  class power_cosine_direction_z : public ray::i_direction<numeric_type> {

  public:
    power_cosine_direction_z(numeric_type exp_) :
      exp(exp_),
      ee(((numeric_type) 2) / (exp + 1)) {}

    util::triple<numeric_type>
    get(rng::i_rng& pRng,
        rng::i_rng::i_state& pRngState1,
        rng::i_rng::i_state& pRngState2
        ) const override final
    {
      auto r1 = ((numeric_type) pRng.get(pRngState)) / ((numeric_type) pRng.max());
      auto r2 = ((numeric_type) pRng.get(pRngState)) / ((numeric_type) pRng.max());
      auto tt = pow(r2, ee);
      auto zz = - sqrtf(tt);
      auto xx = cosf(two_pi * r1) * sqrtf(1 - tt);
      auto yy = sinf(two_pi * r1) * sqrtf(1 - tt);
      auto result = util::triple<numeric_type> {xx, yy, zz};
      util::normalize(result);
      return result;
    }

  private:
    const numeric_type exp;
    const numeric_type ee;
    constexpr static numeric_type two_pi = util::pi() * 2;
  };
}}
