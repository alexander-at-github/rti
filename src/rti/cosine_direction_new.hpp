#pragma once

#include "rti/cos_hemi.hpp"
#include "rti/i_direction.hpp"

namespace rti {
  template<typename T>
  class cosine_direction_new : public i_direction<T> {
  public:
    cosine_direction_new(const rti::triple<rti::triple<T> > pBasis,
                         rti::i_rng* pRng,
                         rti::i_rng::i_state* pRngSeed) :
      mBasis(pBasis),
      mRng(pRng),
      mRngSeed(pRngSeed) {}

    // With unique_ptr one cannot return covariant types
    std::unique_ptr<i_direction<T> > clone() const override final {
      return std::make_unique<cosine_direction_new<T> >(mBasis, mRng, mRngSeed);
    }

    rti::triple<T> get() override final {
      return rti::cos_hemi::get(mBasis, mRng, mRngSeed);
    }

  private:
    rti::triple<rti::triple<T> > mBasis;
    rti::i_rng* mRng;
    rti::i_rng::i_state* mRngSeed;
  };
} // namespace rti
