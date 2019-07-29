#pragma once

#include "rti/ray/cos_hemi.hpp"
#include "rti/ray/i_direction.hpp"

namespace rti { namespace ray {
  template<typename T>
  class cosine_direction : public rti::ray::i_direction<T> {
  public:
    cosine_direction(const rti::util::triple<rti::util::triple<T> > pBasis) :
      mBasis(pBasis) {}

    rti::util::triple<T> get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final {
      return rti::ray::cos_hemi::get(mBasis, pRng, pRngState);
    }

  private:
    rti::util::triple<rti::util::triple<T> > mBasis;
  };
}} // namespace
