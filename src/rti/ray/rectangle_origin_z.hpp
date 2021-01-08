#pragma once

#include "i_origin.hpp"
#include "../util/utils.hpp"

namespace rti { namespace ray {

  template<typename Ty>
  class rectangle_origin_z : public rti::ray::i_origin<Ty> {
  public:

    rectangle_origin_z(Ty pZval, rti::util::pair<Ty> pC1, rti::util::pair<Ty> pC2) :
      mZval(pZval),
      mC1(pC1),
      mC2(pC2) {
      // Rearrange corners if needed such that each value in the pair mC1
      // is smaller or equal to the corresponding value in mC2.
      if (mC1[0] > mC2[0]) rti::util::swap(mC1[0], mC2[0]);
      if (mC1[1] > mC2[1]) rti::util::swap(mC1[1], mC2[1]);
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Condition on ordering of corner points");
    }

    rti::util::triple<Ty> get(rti::rng::i_rng& pRng,
                              rti::rng::i_rng::i_state& pRngState1,
                              rti::rng::i_rng::i_state& pRngState2
                              ) const override final {
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Class invariant on ordering of corner points");
      auto r1 = (Ty) pRng.get(pRngState1);
      auto r2 = (Ty) pRng.get(pRngState2);
      auto rngmin = (Ty) pRng.min();
      auto rngmax = (Ty) pRng.max();
      auto r1prct = (r1 - rngmin) / (rngmax - rngmin + 1); // +1 in order to exclude upper limit !
      auto r2prct = (r2 - rngmin) / (rngmax - rngmin + 1); // +1 in order to exclude upper limit !
      Ty xx = mC1[0] + (mC2[0] - mC1[0]) * r1prct;
      Ty yy = mC1[1] + (mC2[1] - mC1[1]) * r2prct;
      return {xx, yy, mZval};
    }

  private:
    Ty mZval;
    rti::util::pair<Ty> mC1;
    rti::util::pair<Ty> mC2;
  };
}}
