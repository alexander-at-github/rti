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

    rti::util::triple<Ty> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) override final {
      assert(mC1[0] <= mC2[0] && mC1[1] <= mC2[1] && "Class invariant on ordering of corner points");
      auto r1 = (Ty) pRng.get(pRngState);
      auto r2 = (Ty) pRng.get(pRngState);
      auto rngmin = (Ty) pRng.min();
      auto rngmax = (Ty) pRng.max();
      // std::cerr << "r1: " << rngmin << " " << r1 << " " << rngmax << std::endl;
      // std::cerr << "r2: " << rngmin << " " << r2 << " " << rngmax << std::endl;
      auto r1prct = (r1 - rngmin) / (rngmax - rngmin);
      auto r2prct = (r2 - rngmin) / (rngmax - rngmin);
      // std::cerr << r1prct << " " << r2prct << std::endl;
      Ty xx = mC1[0] + (mC2[0] - mC1[0]) * r1prct;
      Ty yy = mC1[1] + (mC2[1] - mC1[1]) * r2prct;
      // std::cerr << mC1[0] << " " << mC2[0] << " " << xx << std::endl;
      // std::cerr << mC1[1] << " " << mC2[1] << " " << yy << std::endl;
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
