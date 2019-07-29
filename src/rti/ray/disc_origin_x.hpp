#pragma once

#include "rti/ray/i_origin.hpp"

namespace rti { namespace ray {
  // T is intended to be a numeric type
  template<typename T>
  class disc_origin_x : public rti::ray::i_origin<T> {
  public:

    disc_origin_x(T pX, T pY, T pZ, T pRadius) :
      mX(pX),
      mY(pY),
      mZ(pZ),
      mRadius(pRadius) {}

    rti::util::triple<T>get(rti::rng::i_rng& pRng, rti::rng::i_rng::i_state& pRngState) const override final {
      T r1 = 1;
      T r2 = 1;
      do {
        r1 = (T) pRng.get(pRngState);
        r2 = (T) pRng.get(pRngState);
        r1 -= pRng.max()/2;
        r2 -= pRng.max()/2;
        r1 = r1 / (pRng.max()/2) * mRadius;
        r2 = r2 / (pRng.max()/2) * mRadius;
        assert(-mRadius <= r1 && r1 <= mRadius &&
               "Error in computin random number in the interval [-mRadius, +mRadius]");
        assert(-mRadius <= r2 && r2 <= mRadius &&
               "Error in computin random number in the interval [-mRadius, +mRadius]");
      } while ( sqrt(r1*r1 + r2*r2) > mRadius );
      return {mX, mY+r1, mZ+r2};
    }

  private:
    T mX;
    T mY;
    T mZ;
    T mRadius;
  };
}} // namespace
