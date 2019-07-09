#pragma once

#include "rti/cos_hemi.hpp"
#include "rti/i_direction.hpp"

namespace rti {
  template<typename T>
  class cosine_direction_new : public i_direction<T> {
  public:
    cosine_direction_new(const rti::triple<rti::triple<T> > pBasis) :
      mBasis(pBasis) {}

    rti::triple<T> get(rti::i_rng& pRng, rti::i_rng::i_state& pRngState) const override final {
      return rti::cos_hemi::get(mBasis, pRng, pRngState);
    }

  private:
    rti::triple<rti::triple<T> > mBasis;
  };
} // namespace rti
