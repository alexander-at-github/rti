#pragma once

#include "rti/ray/i_origin.hpp"

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
      if (mC1[0] > mC2[0]) this->swap(mC1[0], mC2[0]);
      if (mC1[1] > mC2[1]) this->swap(mC1[1], mC2[1]);
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Condition on ordering of corner points");
    }

    rti::util::triple<Ty> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final {
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Class invariant on ordering of corner points");
      auto r1 = (Ty) pRng.get(pRngState);
      auto r2 = (Ty) pRng.get(pRngState);
      auto rngmin = (Ty) pRng.min();
      auto rngmax = (Ty) pRng.max();
      auto r1prct = (rngmax - rngmin) / (r1 - rngmin);
      auto r2prct = (rngmax - rngmin) / (r2 - rngmin);
      Ty xx = mC1[0] + (mC2[0] - mC1[0]) * r1prct;
      Ty yy = mC2[1] + (mC2[1] - mC2[1]) * r2prct;
      return {xx, yy, mZval};
    }

  private:
    Ty mZval;
    rti::util::pair<Ty> mC1;
    rti::util::pair<Ty> mC2;

    void swap(Ty& e1, Ty& e2) const {
      Ty tmp = e2;
      e2 = e1;
      e1 = tmp;
    }
  };
}} // namespace
